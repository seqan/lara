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
#include <set>
#include <tuple>
#include <unordered_map>
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
    std::unordered_map<size_t, size_t> edgeMatching;
    std::vector<PosPair> lines;

    seqan::Rna5String sequenceA;
    seqan::Rna5String sequenceB;

    // the best structural alignment score found so far
    ScoreType bestStructuralAlignmentScore;

    struct EdgeManager
    {
        std::vector<bool> active{};
        size_t size{};
        size_t dim{};

        size_t index(size_t source, size_t target) const
        {
            return dim * source + target;
        }

        size_t source(size_t edgeId) const
        {
            return edgeId / dim;
        }

        size_t target(size_t edgeId) const
        {
            return edgeId % dim;
        }

        bool nonCrossing (size_t idA, size_t idB) const
        {
            bool const smallerA = source(idA) < source(idB) && target(idA) < target(idB);
            bool const smallerB = source(idB) < source(idA) && target(idB) < target(idA);
            return smallerA || smallerB;
        }
    } edges;

    struct InteractionInfo
    {
        ScoreType score;
        size_t dualIdx;
        PriorityQueue::iterator queuePtr;
    };

    std::vector<std::unordered_map<size_t, InteractionInfo>> interaction;

    // every alignment edge holds a priority queue that handles the possible partner edges
    // - the second argument denotes the index of the alignment edges
    // - the first argument holds the actual profit between this pair of alignment edges
    std::vector<PriorityQueue> priorityQ;

    // mapping from the index of the dual variable to the pair of alignment edge indices
    std::vector<PosPair> dualToPairedEdges; // former _YToIndex

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

    void adaptPriorityQ(PosPair pair, ScoreType value)
    {
        priorityQ[pair.first].erase(interaction[pair.first][pair.second].queuePtr);
        interaction[pair.first][pair.second].queuePtr = priorityQ[pair.first].emplace(-value, pair.second).first;
    }

    // return gap score
    ScoreType evaluateLines(std::pair<GappedSeq, GappedSeq> const & alignment, ScoreType gap_open, ScoreType gap_extend)
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
        ScoreType gapScore = 0;

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

    inline ScoreType getSeqScore(SeqScoreMatrix const & mat, size_t idx)
    {
        return sequenceScaleFactor * seqan::score(mat, sequenceA[edges.source(idx)], sequenceB[edges.target(idx)]);
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
        bestStructuralAlignmentScore = -infinity;

        // the following things have to be done in the constructor
        // - given the two RNA structures, determine possible partner edges
        // - compute a set of alignment edges between the RNA structures
        // - initialize the priority queues according to the scores given
        // - provide a mapping between indices (indexing the dual variables
        //   and the actual pair of alignment edges

        edges.active.resize(seqLen.first * seqLen.second);
        edges.size = edges.active.size();
        edges.dim = seqLen.second;
        float const avSeqId = generateEdges(edges.active, sequenceA, sequenceB, params.rnaScore,
                                            static_cast<ScoreType>(params.suboptimalDiff * factor2int));
        sequenceScaleFactor = params.balance * avSeqId + params.sequenceScale;

        priorityQ.resize(edges.size);
        interaction.resize(edges.size);
        dualToPairedEdges.reserve(edges.size);

        // start

        dimension = 0ul;
        for (size_t edgeIdx = 0ul; edgeIdx < edges.size; ++edgeIdx)
        {
            if (!edges.active[edgeIdx])
                continue;

            ScoreType alignScore = getSeqScore(params.rnaScore, edgeIdx);
            priorityQ[edgeIdx].emplace(-alignScore, edgeIdx);

            std::vector<Contact> headContact;
            std::vector<Contact> tailContact;
            extractContacts(headContact, seqan::front(recordA.bppMatrGraphs), edges.source(edgeIdx));
            extractContacts(tailContact, seqan::front(recordB.bppMatrGraphs), edges.target(edgeIdx));

            for (Contact & head : headContact)
            {
                for (Contact & tail : tailContact)
                {
                    size_t partnerIdx = edges.index(head.second, tail.second);
                    if (edges.active[partnerIdx] && edges.nonCrossing(edgeIdx, partnerIdx))
                    {
                        _LOG(3, "     dual idx " << dimension << " = (" << edges.source(edgeIdx)
                                                 << "-" << edges.target(edgeIdx) << ") -> (" << edges.source(partnerIdx)
                                                 << "-" << edges.target(partnerIdx) << ")" << std::endl);

                        ScoreType const structScore = 0.5f * (head.first + tail.first) * factor2int;
                        interaction[edgeIdx][partnerIdx] =
                        {
                            structScore,                                                              // score
                            dimension++,                                                              // dualIdx
                            priorityQ[edgeIdx].emplace(-(structScore + alignScore), partnerIdx).first // queuePtr
                        };
                        dualToPairedEdges.emplace_back(edgeIdx, partnerIdx);
                    }
                }
            }
        }

        // filling the matrix, we're updating the values afterwards, otherwise
        // we had to evaluate _maxProfitScores after every update of either the
        // l or m edge
        for (size_t edgeIdx = 0ul; edgeIdx < edges.size; ++edgeIdx)
        {
            if (edges.active[edgeIdx])
                pssm->set(seqIdx, edges.source(edgeIdx), edges.target(edgeIdx), -priorityQ[edgeIdx].begin()->first);
        }
    }

    void updateScores(std::vector<ScoreType> & dual, std::list<size_t> const & dualIndices, SeqScoreMatrix const & mat)
    {
        for (size_t dualIdx : dualIndices)
        {
            PosPair pair = dualToPairedEdges[dualIdx]; // (l,m)
            ScoreType const newScore = getSeqScore(mat, pair.first) + interaction[pair.first][pair.second].score
                + dual[dualIdx];
            adaptPriorityQ(pair, newScore);
            pssm->set(seqIdx, edges.source(pair.first), edges.target(pair.first),
                      -priorityQ[pair.first].begin()->first);
        }
    }

    ScoreType valid_solution(std::vector<float> & subgradient, std::list<size_t> & subgradientIndices,
                             std::pair<GappedSeq, GappedSeq> const & alignment, unsigned lookahead,
                             SeqScoreMatrix const & mat)
    {
        ScoreType gapScore = evaluateLines(alignment, mat.data_gap_open, mat.data_gap_extend);

        std::vector<size_t> currentStructuralAlignment;
        std::vector<bool> inSolution;
        inSolution.resize(edges.size, false);

        for (PosPair line : lines)
        {
            size_t edgeIdx = edges.index(line.first, line.second);
            SEQAN_ASSERT_MSG(edges.active[edgeIdx], "Alignment match where no alignment edge is defined!");
            currentStructuralAlignment.push_back(edgeIdx);
            inSolution[edgeIdx] = true;
        }

        subgradientIndices.clear();
        for (size_t idx : currentStructuralAlignment)
        {
            // observe all elements with the highest score
            auto fwdIt = priorityQ[idx].begin();
            for (ScoreType const fwdScore = fwdIt->first; fwdScore == fwdIt->first; ++fwdIt)
            {
                // check if there is an interaction cycle: maxProfit(a) == b && maxProfit(b) == a
                bool foundCycle = false;
                size_t const maxPE = fwdIt->second;
                if (inSolution[maxPE])
                {
                    auto revIt = priorityQ[maxPE].begin();
                    for (ScoreType const revScore = revIt->first; !foundCycle && revScore == revIt->first; ++revIt)
                        foundCycle = (revIt->second == idx);
                }

                if (!foundCycle)
                {
                    // adapt dual variables
                    auto dualIt = interaction[idx].find(maxPE);
                    SEQAN_ASSERT_MSG(dualIt != interaction[idx].end(), "Strange condition: undefined interaction.");
                    subgradient[dualIt->second.dualIdx] = 1.0f;
                    subgradientIndices.push_back(dualIt->second.dualIdx);

                    dualIt = interaction[maxPE].find(idx);
                    SEQAN_ASSERT_MSG(dualIt != interaction[maxPE].end(), "Strange condition: undefined interaction.");
                    subgradient[dualIt->second.dualIdx] = -1.0f;
                    subgradientIndices.push_back(dualIt->second.dualIdx);
                }
            }
        }

        ScoreType lowerBound = 0;
        for (size_t idx : currentStructuralAlignment)
            lowerBound += getSeqScore(mat, idx);

        std::unordered_map<size_t, size_t> contacts{};
        if (!subgradientIndices.empty())
        {
            std::vector<std::vector<Contact>> partners{};
            partners.resize(currentStructuralAlignment.size());
            for (size_t idx = 0ul; idx < currentStructuralAlignment.size(); ++idx)
            {
                size_t line = currentStructuralAlignment[idx];
                for (auto const & it : priorityQ[line])
                    if (inSolution[it.second] && line < it.second)
                        partners[idx].emplace_back(interaction[line][it.second].score, it.second);
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
                    lowerBound += interaction[idx][maxPE].score;
                    contacts[idx] = maxPE;
                    contacts[maxPE] = idx;
                }
            }
        }

        // we have to substract the gapcosts, otherwise the lower bound might be higher than the upper bound
        ScoreType primalValue = lowerBound + gapScore;
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
    WeightedAlignedColumns getStructureLines(Parameters const & params, PosPair const & seqIndices) const
    {
        bool swapIdx = seqIndices.first > seqIndices.second;
        WeightedAlignedColumns structureLines{};
        structureLines.first = swapIdx ? PosPair{seqIndices.second, seqIndices.first} : seqIndices;

        if (params.libraryScoreIsLinear)
        {
            auto const mm = std::minmax_element(priorityQ.begin(), priorityQ.end(),
                [] (PriorityQueue const & a, PriorityQueue const & b) { return a.begin()->first > b.begin()->first; });
            ScoreType const minScore = -(mm.first)->begin()->first;
            ScoreType const maxScore = -(mm.second)->begin()->first;
            float const div = 1.f * (params.libraryScoreMax - params.libraryScoreMin) / (maxScore - minScore);

            if (swapIdx)
            {
                for (size_t idx : bestStructuralAlignment)
                {
                    ScoreType val = edgeMatching.count(idx) * (-priorityQ[idx].begin()->first - minScore) * div;
                    structureLines.second.emplace_back(edges.target(idx),
                                                       edges.source(idx),
                                                       params.libraryScoreMin + val);
                }
            }
            else
            {
                for (size_t idx : bestStructuralAlignment)
                {
                    ScoreType val = edgeMatching.count(idx) * (-priorityQ[idx].begin()->first - minScore) * div;
                    structureLines.second.emplace_back(edges.source(idx),
                                                       edges.target(idx),
                                                       params.libraryScoreMin + val);
                }
            }
        }
        else if (swapIdx)
        {
            for (size_t idx : bestStructuralAlignment)
                structureLines.second.emplace_back(edges.target(idx),
                                                   edges.source(idx),
                                                   edgeMatching.count(idx) * 500u + 500u);
        }
        else
        {
            for (size_t idx : bestStructuralAlignment)
                structureLines.second.emplace_back(edges.source(idx),
                                                   edges.target(idx),
                                                   edgeMatching.count(idx) * 500u + 500u);
        }
        return structureLines;
    }
};

} // namespace lara
