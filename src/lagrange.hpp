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
    Parameters params;

    // number of primal/dual variables
    PosPair dimension;

    // mapping of node index to layer
    std::vector<size_t> layer;

    // mapping of node index to position in specific layer
    std::vector<size_t> offset;

    // inverse mapping of layer and offset to  node index
    std::map<PosPair, size_t> invLayerOffset;

    // the arrays hold the alignment scores that are passed onto
    // the pairwiseAlignmentAlgorithm object
    std::vector<std::vector<float>> maxProfitScores;

    // vector holding the maximum profit for each alignment edge
    std::vector<float> maxProfit;

    // vector holding the index of the maximum profit line
    std::vector<size_t> maxProfitEdge;

    std::vector<size_t> bestStructuralAlignment;
    std::map<size_t, size_t> edgeMatching;
    std::vector<PosPair> lines;

    // number of iterations performed so far
    size_t numIterations;

    // number of residues in total
    size_t residueCount;

    seqan::Rna5String sequenceA;
    seqan::Rna5String sequenceB;

    // the best structural alignment score found so far
    float bestStructuralAlignmentScore;
    Alignment bestAlignment;
    Alignment currentAlignment;

    std::vector<size_t> sourceNode;
    std::vector<size_t> targetNode;
    std::map<PosPair, size_t> getEdgeIdx;

    std::map<PosPair, PriorityQueue::iterator> edgeToPriorityQ;

    // primal information used in the branching process
    std::vector<float> primal;

    std::pair<seqan::RnaStructureGraph, seqan::RnaStructureGraph> bppGraphs;
    std::vector<std::vector<Contact>> possiblePartners;

    // scores
    std::vector<float> sequencesScore;
    std::map<PosPair, float> structureScore;

    // every alignment edge holds a priority queue that handles the possible partner edges
    // - the second argument denotes the index of the alignment edges
    // - the first argument holds the actual profit between this pair of alignment edges
    std::vector<PriorityQueue> priorityQ;

    // mapping from pair of edge indices to corresponding dual variable
    std::map<PosPair, size_t> pairedEdgesToDual; // former _IndexToY

    // mapping from the index of the dual variable to the pair of alignment edge indices
    std::vector<PosPair> dualToPairedEdges; // former _YToIndex

    void extractContacts(std::vector<Contact> & contacts, seqan::RnaStructureGraph const & graph, size_t origin)
    {
        for (seqan::RnaAdjacencyIterator adjIt(graph.inter, origin); !seqan::atEnd(adjIt); seqan::goNext(adjIt))
        {
            size_t partner = seqan::value(adjIt);
            float probability = seqan::cargo(seqan::findEdge(graph.inter, origin, partner));
            contacts.emplace_back(partner, probability);
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

    // return gap cost
    float evaluateLines(AlignmentRow const & rowA, AlignmentRow const & rowB)
    {
        typedef typename seqan::Iterator<seqan::Gaps<seqan::Rna5String, seqan::ArrayGaps> const>::Type GapsIterator;

        // Get iterators.
        GapsIterator it0    = seqan::begin(rowA);
        GapsIterator itEnd0 = seqan::end(rowA);
        GapsIterator it1    = seqan::begin(rowB);
        GapsIterator itEnd1 = seqan::end(rowB);

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
                    gapScore += params.laraGapOpen;
                }
                else
                {
                    gapScore += params.laraGapExtend;
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
                    gapScore += params.laraGapOpen;
                }
                else
                {
                    gapScore += params.laraGapExtend;
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

public:
    Lagrange(seqan::RnaRecord const & recordA, seqan::RnaRecord const & recordB, Parameters const & _params)
        : params(_params)
    {
        _LOG(2, recordA.sequence << std::endl << recordB.sequence << std::endl);
        sequenceA = seqan::Rna5String{recordA.sequence};
        sequenceB = seqan::Rna5String{recordB.sequence};
        PosPair seqLen{seqan::length(sequenceA), seqan::length(sequenceB)};

        layer.resize(seqLen.first + seqLen.second);
        offset.resize(seqLen.first + seqLen.second);
        maxProfitScores.resize(seqLen.first);
        for (std::vector<float> & elem : maxProfitScores)
            elem.resize(seqLen.second, negInfinity);

        // initialize()
        numIterations = 0ul;
        residueCount = 0ul;
        dimension.first = 0ul;
        dimension.second = 0ul;

        // score of best alignment
        bestStructuralAlignmentScore = negInfinity;

        // the following things have to be done in the constructor
        // - given the two RNA structures, determine possible partner edges
        // - compute a set of alignment edges between the RNA structures
        // - initialize the priority queues according to the scores given
        // - provide a mapping between indices (indexing the dual variables
        //   and the actual pair of alignment edges

        // create node indices for the RNA objects
        for (size_t idx = 0ul; idx < seqLen.first; ++idx)
        {
            layer[residueCount]                      = 0ul;
            offset[residueCount]                     = idx;
            invLayerOffset[std::make_pair(0ul, idx)] = residueCount++;
        }
        for (size_t idx = 0ul; idx < seqLen.second; ++idx)
        {
            layer[residueCount]                      = 1ul;
            offset[residueCount]                     = idx;
            invLayerOffset[std::make_pair(1ul, idx)] = residueCount++;
        }
        _LOG(2, "residueCount = " << residueCount << " (" << seqLen.first << "+" << seqLen.second << ")" << std::endl);

        generateEdges(getEdgeIdx, sequenceA, sequenceB, params.laraScoreMatrix, params.suboptimalDiff);
        size_t numEdges = getEdgeIdx.size();

        sourceNode.reserve(numEdges);
        targetNode.reserve(numEdges);
        primal.resize(numEdges, 0.0);
        sequencesScore.resize(numEdges, 0.0);
        possiblePartners.resize(numEdges);
        maxProfit.resize(numEdges, negInfinity);
        maxProfitEdge.resize(numEdges, SIZE_MAX);
        priorityQ.resize(numEdges);
        dualToPairedEdges.reserve(numEdges);

        for (std::pair<PosPair const, size_t> const & edge : getEdgeIdx)
        {
            sourceNode.push_back(edge.first.first);
            targetNode.push_back(edge.first.second);
        }
        bppGraphs = std::make_pair(seqan::front(recordA.bppMatrGraphs), seqan::front(recordB.bppMatrGraphs));
        start();
    }

    void start()
    {
        size_t dualIdx = 0ul;
        size_t numPartners = 0ul;
        for (size_t edgeIdx = 0ul; edgeIdx < sourceNode.size(); ++edgeIdx)
        {
            size_t headNode = sourceNode[edgeIdx];
            size_t tailNode = targetNode[edgeIdx];
            std::vector<Contact> headContact;
            std::vector<Contact> tailContact;
            extractContacts(headContact, bppGraphs.first, headNode);
            extractContacts(tailContact, bppGraphs.second, tailNode);

            float alignScore = seqan::score(params.laraScoreMatrix,
                                            sequenceA[sourceNode[edgeIdx]],
                                            sequenceB[targetNode[edgeIdx]]);

            possiblePartners[edgeIdx].emplace_back(edgeIdx, alignScore);
            structureScore[std::make_pair(edgeIdx, edgeIdx)] = negInfinity;
            sequencesScore[edgeIdx] = alignScore;
            priorityQ[edgeIdx].emplace(-alignScore, edgeIdx);
            maxProfit[edgeIdx] = alignScore;
            maxProfitEdge[edgeIdx] = edgeIdx;

            if (!headContact.empty() && !tailContact.empty()) // there are interactions
            {
                for (Contact & head : headContact)
                {
                    for (Contact & tail : tailContact)
                    {
                        auto partnerIter = getEdgeIdx.find(std::make_pair(head.first, tail.first));
                        if (partnerIter != getEdgeIdx.end() && nonCrossingEdges(edgeIdx, partnerIter->second))
                        {
                            size_t partnerIdx = partnerIter->second;
                            PosPair interaction{edgeIdx, partnerIdx};
                            float structScore = 0.5f * (head.second + tail.second);
                            possiblePartners[edgeIdx].emplace_back(partnerIdx, structScore);
                            structureScore[interaction] = structScore;
                            sequencesScore[edgeIdx] = alignScore;

                            // insert element into the priority queue
                            auto res = priorityQ[edgeIdx].emplace(-(structScore + alignScore), partnerIdx);
                            edgeToPriorityQ[interaction] = res.first;
                            _LOG(2, "dual idx " << dualIdx << " = (" << sourceNode[edgeIdx]
                                                    << "-" << targetNode[edgeIdx] << ") -> (" << sourceNode[partnerIdx]
                                                    << "-" << targetNode[partnerIdx] << ")" << std::endl);
                            pairedEdgesToDual[interaction] = dualIdx++;
                            dualToPairedEdges.emplace_back(edgeIdx, partnerIdx);

                            if (structScore + alignScore > maxProfit[edgeIdx])
                            {
                                _LOG(2, "maxProfit[" << edgeIdx << " (" << head.first << "," << tail.first
                                                         << ")] = " << structScore << " + " << alignScore
                                                         << " \tpartner " << partnerIdx << std::endl);
                                maxProfit[edgeIdx] = structScore + alignScore;
                                maxProfitEdge[edgeIdx] = partnerIdx;
                            }
                        }
                    }
                }
            }
            numPartners += possiblePartners[edgeIdx].size();
        }
        _LOG(2, "Average number of partner edges = " << 1.0 * numPartners / sourceNode.size() << std::endl);
        dimension.second = dualIdx; // number of interactions that are observed
        dimension.first = sourceNode.size(); // number of alignment edges (lines)
    }

    void updateScores(std::vector<float> & dual, std::list<size_t> const & dualIndices)
    {
        for (size_t dualIdx : dualIndices)
        {
            PosPair pair = dualToPairedEdges[dualIdx]; // (l,m)
            adaptPriorityQ(pair, sequencesScore[pair.first] + structureScore[pair] + dual[dualIdx]);
            assert(!priorityQ.empty());
            auto maxElement = priorityQ[pair.first].begin();
            maxProfit[pair.first] = -maxElement->first;     // negative priority value (maximum weight)
            maxProfitEdge[pair.first] = maxElement->second; // information
        }

        // filling the matrix, we're updating the values afterwards, otherwise
        // we had to evaluate _maxProfitScores after every update of either the
        // l or m edge
        for (size_t edgeIdx = 0ul; edgeIdx < sourceNode.size(); ++edgeIdx)
        {
            maxProfitScores[sourceNode[edgeIdx]][targetNode[edgeIdx]] = maxProfit[edgeIdx];
        }

        _LOG(3, "maxProfitScores" << std::endl);
        for (auto & row : maxProfitScores)
        {
            _LOG(3, "[ ");
            for (float sc : row)
            {
                _LOG(3, std::setw(14) << sc << " ");
            }
            _LOG(3, "]" << std::endl);
        }
    }

    /*!
     * \brief Performs the structural alignment.
     * \return The dual value (upper bound, solution of relaxed problem).
     */
    float relaxed_solution(float gapOpen, float gapExtend)
    {
        seqan::Score<float, seqan::RnaStructureScore> scoreAdaptor(&maxProfitScores, gapOpen, gapExtend);
        seqan::resize(seqan::rows(currentAlignment), 2);
        seqan::assignSource(seqan::row(currentAlignment, 0), sequenceA);
        seqan::assignSource(seqan::row(currentAlignment, 1), sequenceB);

        // perform the alignment
        return seqan::globalAlignment(currentAlignment, scoreAdaptor, seqan::AffineGaps());
    }

    float valid_solution(std::vector<float> & subgradient, std::list<size_t> & subgradientIndices)
    {
        float gapScore = evaluateLines(seqan::row(currentAlignment, 0), seqan::row(currentAlignment, 1));

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
            ++primal[edgeIdx]; // primal counter
        }

        subgradientIndices.clear();
        for (size_t idx : currentStructuralAlignment)
        {
            size_t const & maxPE = maxProfitEdge[idx];

            if (inSolution[maxPE] && maxProfitEdge[maxPE] == idx)
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
        std::map<size_t, size_t> contacts{};
        if (!subgradientIndices.empty())
        {
            Matching mwm(sequencesScore, possiblePartners, params.matching);
            lowerBound = mwm.computeScore(currentStructuralAlignment, inSolution);
            contacts = mwm.getContacts();
        }
        else
        {
            for (size_t idx : currentStructuralAlignment)
            {
                size_t const & maxPE = maxProfitEdge[idx];
                lowerBound += sequencesScore[idx];
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
            size_t const & maxPE = maxProfitEdge[idx];
            _LOG(3, "Alignment[" << idx << "; " << sourceNode[idx] << "," << targetNode[idx] << "] maxProfitEdge ["
                                      << maxPE << "; " << sourceNode[maxPE] << "," << targetNode[maxPE] << "] score "
                                      << maxProfit[idx] << " inSolution " << inSolution[maxPE] << " rec "
                                      << (maxProfitEdge[maxPE] == idx) << std::endl);
        }

        // we have to substract the gapcosts, otherwise the lower bound might be higher than the upper bound
        float primalValue = lowerBound + gapScore;
        _LOG(2, "primal " << primalValue << " = " << lowerBound << " (lb) + " << gapScore << " (gap)" << std::endl);

        // store the best alignment found so far
        if (primalValue > bestStructuralAlignmentScore)
        {
            bestStructuralAlignmentScore = primalValue;
            bestStructuralAlignment = currentStructuralAlignment;
            edgeMatching = contacts;
            bestAlignment = Alignment(currentAlignment);
        }
        ++numIterations;
        return primalValue;
    }

    PosPair getDimension()
    {
        return dimension;
    }

    Alignment & getAlignment()
    {
        return bestAlignment;
    }

    std::vector<std::tuple<size_t, size_t, bool>> getStructureLines() const
    {
        std::vector<std::tuple<size_t, size_t, bool>> structureLines{};
        for (size_t idx : bestStructuralAlignment)
        {
            structureLines.emplace_back(sourceNode[idx] + 1, targetNode[idx] + 1, (edgeMatching.count(idx) == 1));
        }
        return structureLines;
    }
};

} // namespace lara

