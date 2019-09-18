// ===========================================================================
//                LaRA: Lagrangian Relaxed structural Alignment
// ===========================================================================
// Copyright (c) 2016-2019, Jörg Winkler, Freie Universität Berlin
// Copyright (c) 2016-2019, Gianvito Urgese, Politecnico di Torino
// Copyright (c) 2006-2019, Knut Reinert, Freie Universität Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
// * Neither the name of Jörg Winkler, Gianvito Urgese, Knut Reinert,
//   the FU Berlin or the Politecnico di Torino nor the names of
//   its contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#pragma once

/*!\file lagrange.hpp
 * \brief This file contains data structures and algorithms for lagrange relaxation.
 */

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include <seqan/rna_io.h>
#include <seqan/score.h>

#include "data_types.hpp"
#include "edge_filter.hpp"
#include "parameters.hpp"
#include "score.hpp"
#include "matching.hpp"

namespace lara
{

class Lagrange
{
private:
    // number of dual variables = number of interactions that are observed
    size_t dimension;

    std::vector<size_t> bestStructuralAlignment;
    std::map<size_t, size_t> edgeMatching;
    std::vector<PosPair> lines;

    seqan::Rna5String sequenceA;
    seqan::Rna5String sequenceB;

    // the best structural alignment score found so far
    float bestStructuralAlignmentScore;
    //Alignment bestAlignment;

    std::vector<size_t> sourceNode;
    std::vector<size_t> targetNode;
    std::map<PosPair, size_t> getEdgeIdx;

    std::pair<seqan::RnaStructureGraph, seqan::RnaStructureGraph> bppGraphs;

    // scores
    std::map<PosPair, float> structureScore;

    // every alignment edge holds a priority queue that handles the possible partner edges
    // - the second argument denotes the index of the alignment edges
    // - the first argument holds the actual profit between this pair of alignment edges
    std::vector<PriorityQueue> priorityQ;
    std::map<PosPair, PriorityQueue::iterator> edgeToPriorityQ;

    // mapping from pair of edge indices to corresponding dual variable
    std::map<PosPair, size_t> pairedEdgesToDual; // former _IndexToY

    // mapping from the index of the dual variable to the pair of alignment edge indices
    std::vector<PosPair> dualToPairedEdges; // former _YToIndex

    std::pair<unsigned, unsigned> libraryScore{};
    bool libraryScoreIsLinear{};

    float gap_open;
    float gap_extend;

    RnaScoreType * pssm;
    size_t seqIdx;
    float sequenceScaleFactor;

    static void extractContacts(std::vector<Contact> & contacts, seqan::RnaStructureGraph const & graph, size_t origin)
    {
        for (seqan::RnaAdjacencyIterator adjIt(graph.inter, origin); !seqan::atEnd(adjIt); seqan::goNext(adjIt))
        {
            size_t partner = seqan::value(adjIt);
            float probability = seqan::cargo(seqan::findEdge(graph.inter, origin, partner));
            contacts.emplace_back(probability, partner);
        }
    }

    bool nonCrossingEdges (size_t edgeA, size_t edgeB)
    {
        bool smallerA = sourceNode[edgeA] < sourceNode[edgeB] && targetNode[edgeA] < targetNode[edgeB];
        bool smallerB = sourceNode[edgeB] < sourceNode[edgeA] && targetNode[edgeB] < targetNode[edgeA];
        return smallerA || smallerB;
    }

    void adaptPriorityQ(PosPair pair, float value)
    {
        priorityQ[pair.first].erase(edgeToPriorityQ[pair]);
        edgeToPriorityQ[pair] = priorityQ[pair.first].emplace(-value, pair.second).first;
    }

    // return gap score
    float evaluateLines(std::pair<GappedSeq, GappedSeq> const & alignment)
    {
        // Get iterators.
        auto it0    = seqan::begin(alignment.first);
        auto itEnd0 = seqan::end(alignment.first);
        auto it1    = seqan::begin(alignment.second);
        auto itEnd1 = seqan::end(alignment.second);

        // State whether we have already opened a gap.
        bool isGapOpen0 = false;
        bool isGapOpen1 = false;

        // Keep track of the sequence positions for lines.
        PosPair sourcePos(0ul, 0ul);
        lines.clear();

        // Sum up gap score.
        float gapScore = 0.0f;

        while (it0 != itEnd0 && it1 != itEnd1)
        {
            // Gaps in first sequence
            if (seqan::isGap(it0))
            {
                if (!isGapOpen0)
                {
                    gapScore += gap_open;
                }
                else
                {
                    gapScore += gap_extend;
                }
                isGapOpen0 = true;
                ++sourcePos.second;
            }
            else
            {
                isGapOpen0 = false;
            }

            // Gaps in second sequence
            if (seqan::isGap(it1))
            {
                if (!isGapOpen1)
                {
                    gapScore += gap_open;
                }
                else
                {
                    gapScore += gap_extend;
                }
                isGapOpen1 = true;
                ++sourcePos.first;
            }
            else
            {
                isGapOpen1 = false;
            }

            // Match or mismatch
            if (!seqan::isGap(it0) && !seqan::isGap(it1))
            {
                // create a line
                lines.push_back(sourcePos);
                ++sourcePos.first;
                ++sourcePos.second;
            }

            ++it0;
            ++it1;
        }
        SEQAN_ASSERT(it0 == itEnd0);
        SEQAN_ASSERT(it1 == itEnd1);
        return gapScore;
    }

    inline float getSeqScore(RnaScoreMatrix const & mat, size_t edgeIdx)
    {
        return sequenceScaleFactor * seqan::score(mat, sequenceA[sourceNode[edgeIdx]], sequenceB[targetNode[edgeIdx]]);
    }

public:
    Lagrange(seqan::RnaRecord const & recordA, seqan::RnaRecord const & recordB,
             Parameters const & params, RnaScoreType * score, size_t sidx) : pssm(score), seqIdx(sidx)
    {
        _LOG(3, "     " << recordA.sequence << "\n     " << recordB.sequence << std::endl);
        sequenceA = seqan::Rna5String{recordA.sequence};
        sequenceB = seqan::Rna5String{recordB.sequence};
        PosPair seqLen{seqan::length(sequenceA), seqan::length(sequenceB)};

        // score of best alignment
        bestStructuralAlignmentScore = negInfinity;

        // the following things have to be done in the constructor
        // - given the two RNA structures, determine possible partner edges
        // - compute a set of alignment edges between the RNA structures
        // - initialize the priority queues according to the scores given
        // - provide a mapping between indices (indexing the dual variables
        //   and the actual pair of alignment edges

        float const avSeqId = generateEdges(getEdgeIdx, sequenceA, sequenceB, params.rnaScore,
                                            params.suboptimalDiff);
        sequenceScaleFactor = params.balance * avSeqId + params.sequenceScale;
        size_t numEdges = getEdgeIdx.size();

        sourceNode.reserve(numEdges);
        targetNode.reserve(numEdges);
        priorityQ.resize(numEdges);
        dualToPairedEdges.reserve(numEdges);

        for (std::pair<PosPair const, size_t> const & edge : getEdgeIdx)
        {
            sourceNode.push_back(edge.first.first);
            targetNode.push_back(edge.first.second);
        }
        bppGraphs = std::make_pair(seqan::front(recordA.bppMatrGraphs), seqan::front(recordB.bppMatrGraphs));
        libraryScore = std::make_pair(params.libraryScoreMin, params.libraryScoreMax);
        libraryScoreIsLinear = params.libraryScoreIsLinear;

        // start

        dimension = 0ul;
        for (size_t edgeIdx = 0ul; edgeIdx < sourceNode.size(); ++edgeIdx)
        {
            size_t headNode = sourceNode[edgeIdx];
            size_t tailNode = targetNode[edgeIdx];
            std::vector<Contact> headContact;
            std::vector<Contact> tailContact;
            extractContacts(headContact, bppGraphs.first, headNode);
            extractContacts(tailContact, bppGraphs.second, tailNode);

            float alignScore = getSeqScore(params.rnaScore, edgeIdx);

            structureScore[std::make_pair(edgeIdx, edgeIdx)] = negInfinity;
            priorityQ[edgeIdx].emplace(-alignScore, edgeIdx);

            if (!headContact.empty() && !tailContact.empty()) // there are interactions
            {
                for (Contact & head : headContact)
                {
                    for (Contact & tail : tailContact)
                    {
                        auto partnerIter = getEdgeIdx.find(std::make_pair(head.second, tail.second));
                        if (partnerIter != getEdgeIdx.end() && nonCrossingEdges(edgeIdx, partnerIter->second))
                        {
                            size_t partnerIdx = partnerIter->second;
                            PosPair interaction{edgeIdx, partnerIdx};
                            float structScore = 0.5f * (head.first + tail.first);
                            structureScore[interaction] = structScore;

                            // insert element into the priority queue
                            auto res = priorityQ[edgeIdx].emplace(-(structScore + alignScore), partnerIdx);
                            edgeToPriorityQ[interaction] = res.first;
                            _LOG(3, "     dual idx " << dimension << " = (" << sourceNode[edgeIdx]
                                                     << "-" << targetNode[edgeIdx] << ") -> (" << sourceNode[partnerIdx]
                                                     << "-" << targetNode[partnerIdx] << ")" << std::endl);
                            pairedEdgesToDual[interaction] = dimension++;
                            dualToPairedEdges.emplace_back(edgeIdx, partnerIdx);
                        }
                    }
                }
            }
        }

        // filling the matrix, we're updating the values afterwards, otherwise
        // we had to evaluate _maxProfitScores after every update of either the
        // l or m edge
        gap_open = params.rnaScore.data_gap_open;
        gap_extend = params.rnaScore.data_gap_extend;
        for (size_t edgeIdx = 0ul; edgeIdx < sourceNode.size(); ++edgeIdx)
        {
            pssm->set(seqIdx, sourceNode[edgeIdx], targetNode[edgeIdx],
                      static_cast<int32_t>(-priorityQ[edgeIdx].begin()->first * factor2int));
        }
    }

    void updateScores(std::vector<float> & dual, std::list<size_t> const & dualIndices, RnaScoreMatrix const & mat)
    {
        for (size_t dualIdx : dualIndices)
        {
            PosPair pair = dualToPairedEdges[dualIdx]; // (l,m)
            adaptPriorityQ(pair, getSeqScore(mat, pair.first) + structureScore[pair] + dual[dualIdx]);
            assert(!priorityQ.empty());
            pssm->set(seqIdx, sourceNode[pair.first], targetNode[pair.first],
                      static_cast<int32_t>(-priorityQ[pair.first].begin()->first * factor2int));
        }
    }

    float valid_solution(std::vector<float> & subgradient, std::list<size_t> & subgradientIndices,
                         std::pair<GappedSeq, GappedSeq> const & alignment, unsigned lookahead,
                         RnaScoreMatrix const & mat)
    {
        float gapScore = evaluateLines(alignment);

        std::vector<size_t> currentStructuralAlignment;
        std::vector<bool> inSolution;
        inSolution.resize(sourceNode.size(), false);

        for (PosPair line : lines)
        {
            auto edgeIdxIt = getEdgeIdx.find(line);
            SEQAN_ASSERT_MSG(edgeIdxIt != getEdgeIdx.end(), "Alignment match where no alignment edge is defined!");

            size_t edgeIdx = edgeIdxIt->second;
            currentStructuralAlignment.push_back(edgeIdx);
            inSolution[edgeIdx] = true;
        }

        subgradientIndices.clear();
        for (size_t idx : currentStructuralAlignment)
        {
            size_t const & maxPE = priorityQ[idx].begin()->second;

            if (inSolution[maxPE] && priorityQ[maxPE].begin()->second == idx)
                continue;

            auto dualIt = pairedEdgesToDual.find(std::make_pair(idx, maxPE));
            if (dualIt != pairedEdgesToDual.end())
            {
                subgradient[dualIt->second] = 1.0f;
                subgradientIndices.push_back(dualIt->second);
            }
            else
            {
                std::cerr << "Strange condition: pairedEdgesToDual(" << idx << "," << maxPE
                          << ") not defined!" << std::endl;
                exit(1);
            }

            dualIt = pairedEdgesToDual.find(std::make_pair(maxPE, idx));
            if (dualIt != pairedEdgesToDual.end())
            {
                subgradient[dualIt->second] = -1.0f;
                subgradientIndices.push_back(dualIt->second);
            }
            else
            {
                std::cerr << "Strange condition: pairedEdgesToDual(" << maxPE << "," << idx
                          << ") not defined!" << std::endl;
                exit(1);
            }
        }

        float lowerBound = 0.0f;
        for (size_t idx : currentStructuralAlignment)
            lowerBound += getSeqScore(mat, idx);

        std::map<size_t, size_t> contacts{};
        if (!subgradientIndices.empty())
        {
            std::vector<std::vector<Contact>> partners{};
            partners.resize(currentStructuralAlignment.size());
            for (size_t idx = 0ul; idx < currentStructuralAlignment.size(); ++idx)
            {
                size_t line = currentStructuralAlignment[idx];
                for (auto const & it : priorityQ[line])
                    if (inSolution[it.second] && line < it.second)
                        partners[idx].emplace_back(structureScore[std::make_pair(line, it.second)], it.second);
            }

            Matching mwm(partners, lookahead);
            lowerBound += mwm.computeScore(currentStructuralAlignment);
            contacts = mwm.getContacts();
        }
        else
        {
            for (size_t idx : currentStructuralAlignment)
            {
                size_t const & maxPE = priorityQ[idx].begin()->second;
                if (idx != maxPE)
                {
                    lowerBound += structureScore[PosPair(idx, maxPE)];
                    contacts[idx] = maxPE;
                    contacts[maxPE] = idx;
                }
            }
        }

        for (size_t idx : currentStructuralAlignment)
        {
            size_t const & maxPE = priorityQ[idx].begin()->second;
            _LOG(3, "     Alignment[" << idx << "; " << sourceNode[idx] << "," << targetNode[idx] << "] maxProfitEdge ["
                                      << maxPE << "; " << sourceNode[maxPE] << "," << targetNode[maxPE] << "] score "
                                      << -priorityQ[idx].begin()->first << " inSolution " << inSolution[maxPE]
                                      << " rec " << (priorityQ[maxPE].begin()->second == idx) << std::endl);
        }

        // we have to substract the gapcosts, otherwise the lower bound might be higher than the upper bound
        float primalValue = lowerBound + gapScore;
        _LOG(3, "     primal " << primalValue << " = " << lowerBound << " (lb) + " << gapScore << " (gp)" << std::endl);

        // store the best alignment found so far
        if (primalValue > bestStructuralAlignmentScore)
        {
            bestStructuralAlignmentScore = primalValue;
            bestStructuralAlignment = currentStructuralAlignment;
            edgeMatching = contacts;
        }
        return primalValue;
    }

    size_t getDimension()
    {
        return dimension;
    }

    /*!
     * \brief Calculate the scores for the T-Coffee library.
     * \return A vector of triples of position, position and score.
     * \details
     * The score is normalised to [l_min..l_max] if the pair is included in the matching
     * and equals l_min otherwise. The pair {l-min, l_max} is the libscore parameter of LaRA.
     */
    std::vector<std::tuple<size_t, size_t, unsigned>> getStructureLines() const
    {
        std::vector<std::tuple<size_t, size_t, unsigned>> structureLines{};
        if (libraryScoreIsLinear)
        {
            auto const mm = std::minmax_element(priorityQ.begin(), priorityQ.end(),
                [] (PriorityQueue const & a, PriorityQueue const & b) { return a.begin()->first > b.begin()->first; });
            float const minScore = -(mm.first)->begin()->first;
            float const maxScore = -(mm.second)->begin()->first;
            float const div = 1. * (libraryScore.second - libraryScore.first) / (maxScore - minScore);
            for (size_t idx : bestStructuralAlignment)
            {
                auto val = static_cast<unsigned>(edgeMatching.count(idx) * (-priorityQ[idx].begin()->first - minScore)
                                                                         * div);
                structureLines.emplace_back(sourceNode[idx] + 1, targetNode[idx] + 1, libraryScore.first + val);
            }
        }
        else
        {
            for (size_t idx : bestStructuralAlignment)
                structureLines.emplace_back(sourceNode[idx] + 1, targetNode[idx] + 1,
                                            edgeMatching.count(idx) * 500u + 500u);
        }
        return structureLines;
    }
};

} // namespace lara
