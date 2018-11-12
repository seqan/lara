// ===========================================================================
//                LaRA: Lagrangian Relaxed structural Alignment
// ===========================================================================
// Copyright (c) 2016-2018, Jörg Winkler, Freie Universität Berlin
// Copyright (c) 2016-2018, Gianvito Urgese, Politecnico di Torino
// Copyright (c) 2006-2018, Knut Reinert, Freie Universität Berlin
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

#include <iostream>
#include <map>
#include <cstdint>
#include <utility>
#include <vector>

#ifdef LEMON_FOUND
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/concepts/graph.h>
#endif

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include <seqan/rna_io.h>
#include <seqan/score.h>

#include "data_types.hpp"
#include "edge_filter.hpp"
#include "parameters.hpp"
#include "score.hpp"

namespace lara
{

class Lagrange
{
private:
    Parameters & params;

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
    std::vector<std::vector<double>> maxProfitScores;

    // vector holding the maximum profit for each alignment edge
    std::vector<double> maxProfit;

    // vector holding the index of the maximum profit line
    std::vector<size_t> maxProfitEdge;

    // holds the initial maxProfitScores as set in StartUpLagrange()
    std::vector<double>              bestLagrangianMultipliers;
    std::vector<std::vector<double>> bestUpperBoundScores;

    std::vector<size_t> bestStructuralAlignment;
    std::map<size_t, size_t> edgeMatching;
    std::vector<PosPair> lines;

    // number of iterations performed so far
    size_t numIterations;

    // number of residues in total
    size_t residueCount;

    // best upper bound found so far
    std::vector<double> allUpperBounds;

    seqan::Rna5String sequenceA;
    seqan::Rna5String sequenceB;

    bool doMatching; // True for bpp matrix input

    // the best structural alignment score found so far
    double bestStructuralAlignmentScore;
    Alignment bestAlignment;

    // the best (lowest) upper bound found so far
    double bestUpperBoundScore;

    std::vector<size_t> sourceNode;
    std::vector<size_t> targetNode;
    std::map<PosPair, size_t> getEdgeIdx;

    std::map<PosPair, PriorityQueue::iterator> edgeToPriorityQ;

    // primal information used in the branching process
    std::vector<double> primal;

    std::pair<seqan::RnaStructureGraph, seqan::RnaStructureGraph> bppGraphs;
    std::vector<std::vector<Contact>> possiblePartners;

    // scores
    std::vector<double> sequencesScore;
    std::map<PosPair, double> structureScore;

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
            double probability = seqan::cargo(seqan::findEdge(graph.inter, origin, partner));
            contacts.emplace_back(partner, probability);
        }
    }

    bool nonCrossingEdges (size_t edgeA, size_t edgeB)
    {
        bool smallerA = sourceNode[edgeA] < sourceNode[edgeB] && targetNode[edgeA] < targetNode[edgeB];
        bool smallerB = sourceNode[edgeB] < sourceNode[edgeA] && targetNode[edgeB] < targetNode[edgeA];
        return smallerA || smallerB;
    }

    void adaptPriorityQ(PosPair pair, double value)
    {
        priorityQ[pair.first].erase(edgeToPriorityQ[pair]);
        edgeToPriorityQ[pair] = priorityQ[pair.first].emplace(-value, pair.second).first;
    }

    // return gap cost
    double evaluateLines(AlignmentRow const & rowA, AlignmentRow const & rowB)
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
        double gapScore = 0.0;

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

#ifdef LEMON_FOUND
    double computeMatching(std::map<size_t, size_t> & contacts,
                           std::vector<size_t> const & currentAlignment,
                           std::vector<bool> const & inSolution)
    {
        double score = 0.0;
        contacts.clear();

        lemon::SmartGraph lemonG;
        std::map<size_t, lemon::SmartGraph::Node> nodes;

        for (size_t const & line : currentAlignment)
        {
            contacts[line] = line;
            nodes[line] = lemonG.addNode();
        }

        typedef lemon::SmartGraph::EdgeMap<double> EdgeMap;
        EdgeMap weight(lemonG);
        lemon::SmartGraph::EdgeMap<PosPair> interactions(lemonG);
        std::map<PosPair, bool> computed;
        for (size_t const & line : currentAlignment)
        {
            score += sequencesScore[line];
            for (Contact const & contact : possiblePartners[line])
            {
                PosPair interaction{line, contact.first};
                auto res = computed.find(interaction);
                if (res == computed.end() && inSolution[contact.first] && line != contact.first)
                {
                    PosPair revInteraction = std::make_pair(contact.first, line);
                    _VVV(params, "addEdge " << line << " " << contact.first << " ("
                                            << structureScore[interaction] + structureScore[revInteraction] << ")");
                    auto newEdge = lemonG.addEdge(nodes[line], nodes[contact.first]);
                    weight[newEdge] = structureScore[interaction] + structureScore[revInteraction];
                    interactions[newEdge] = interaction;
                    computed[interaction] = true;
                    computed[revInteraction] = true;
                }
            }
        }
        lemon::MaxWeightedMatching<lemon::SmartGraph, EdgeMap> mwm(lemonG, weight);
        mwm.run();
        double lemonweight = mwm.matchingWeight();
        for (lemon::SmartGraph::EdgeIt edgeIt(lemonG); edgeIt!=lemon::INVALID; ++edgeIt)
        {
            PosPair inter = interactions[edgeIt];
            if (mwm.matching(edgeIt))
            {
                contacts[inter.first] = inter.second;
                contacts[inter.second] = inter.first;
                _VVV(params, "lemon matches  " << inter.first << " " << inter.second << " " << weight[edgeIt]);
            }
            else
            {
                _VVV(params, "lemon excludes " << inter.first << " " << inter.second << " " << weight[edgeIt]);
            }
        }
        _VV(params, "\nlower bound: seq " << score << " + str " << lemonweight);
        return score + lemonweight;
    }
#endif



public:
    Lagrange(seqan::RnaRecord const & recordA, seqan::RnaRecord const & recordB, Parameters & _params) : params(_params)
    {
        _VV(params, recordA.sequence << std::endl << recordB.sequence);
        sequenceA = seqan::Rna5String{recordA.sequence};
        sequenceB = seqan::Rna5String{recordB.sequence};
        PosPair seqLen{seqan::length(sequenceA), seqan::length(sequenceB)};

        layer.resize(seqLen.first + seqLen.second);
        offset.resize(seqLen.first + seqLen.second);
        maxProfitScores.resize(seqLen.first);
        for (std::vector<double> & elem : maxProfitScores)
            elem.resize(seqLen.second, negInfinity);
        bestLagrangianMultipliers.resize(seqLen.first * seqLen.second);
        bestUpperBoundScores.resize(seqLen.first);
        for (std::vector<double> & elem : bestUpperBoundScores)
            elem.resize(seqLen.second, negInfinity);

        // initialize()
        numIterations = 0ul;
        residueCount = 0ul;
        dimension.first = 0ul;
        dimension.second = 0ul;
        doMatching = false;

        // score of best alignment
        bestStructuralAlignmentScore = negInfinity;
        bestUpperBoundScore = posInfinity;

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
        _VV(params, "residueCount = " << residueCount << " (" << seqLen.first << "+" << seqLen.second << ")");

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

            double alignScore = seqan::score(params.laraScoreMatrix,
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
                            double structScore = 0.5 * (head.second + tail.second);
                            possiblePartners[edgeIdx].emplace_back(partnerIdx, structScore);
                            structureScore[interaction] = structScore;
                            sequencesScore[edgeIdx] = alignScore;

                            // insert element into the priority queue
                            auto res = priorityQ[edgeIdx].emplace(-(structScore + alignScore), partnerIdx);
                            edgeToPriorityQ[interaction] = res.first;
                            _VV(params, "dual idx " << dualIdx << " = (" << sourceNode[edgeIdx]
                                                    << "-" << targetNode[edgeIdx] << ") -> (" << sourceNode[partnerIdx]
                                                    << "-" << targetNode[partnerIdx] << ")");
                            pairedEdgesToDual[interaction] = dualIdx++;
                            dualToPairedEdges.emplace_back(edgeIdx, partnerIdx);

                            if (structScore + alignScore > maxProfit[edgeIdx])
                            {
                                _VV(params, "maxProfit[" << edgeIdx << " (" << head.first << "," << tail.first
                                                         << ")] = " << structScore << " + " << alignScore
                                                         << " \tpartner " << partnerIdx);
                                maxProfit[edgeIdx] = structScore + alignScore;
                                maxProfitEdge[edgeIdx] = partnerIdx;
                            }
                        }
                    }
                }
            }

            doMatching = doMatching || (possiblePartners[edgeIdx].size() > 2);
            numPartners += possiblePartners[edgeIdx].size();
        }
        _VV(params, "Average number of partner edges = " << 1.0 * numPartners / sourceNode.size());
        dimension.second = dualIdx;
        dimension.first = sourceNode.size();
    }

    void evaluate(std::vector<double> & dual,
                  std::list<size_t> & dualIndices,
                  double & dualValue,   // upper bound
                  double & primalValue, // lower bound
                  std::vector<double> & subgradient,
                  std::list<size_t> & subgradientIndices)
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

        if (params.verbose >= 3)
        {
            std::cerr << "maxProfitScores" << std::endl;
            for (auto & row : maxProfitScores)
            {
                std::cerr << "[ ";
                for (double sc : row)
                {
                    std::cerr << std::setw(14) << sc << " ";
                }
                std::cerr << "]" << std::endl;
            }
        }

        seqan::Score<double, seqan::RnaStructureScore> scoreAdaptor(&maxProfitScores,
                                                                    params.laraGapOpen,
                                                                    params.laraGapExtend);
        Alignment alignment;
        seqan::resize(seqan::rows(alignment), 2);
        AlignmentRow & rowA = seqan::row(alignment, 0);
        AlignmentRow & rowB = seqan::row(alignment, 1);
        _VV(params, "sequence length " << length(sequenceA) << " " << length(sequenceB));

        seqan::assignSource(rowA, sequenceA);
        seqan::assignSource(rowB, sequenceB);

        // perform the alignment
        double optScore = seqan::globalAlignment(alignment, scoreAdaptor, seqan::AffineGaps());
        double gapScore = evaluateLines(rowA, rowB);

        std::vector<size_t> currentStructuralAlignment;
        std::vector<bool> inSolution;
        inSolution.resize(sourceNode.size(), false);
        double scoreseq = gapScore;

        for (PosPair line : lines)
        {
            //double sc = seqan::score(scoreAdaptor, line.first, line.second);
            double sc = maxProfitScores[line.first][line.second];
            scoreseq += sc;

            auto edgeIdxIt = getEdgeIdx.find(line);
            if (edgeIdxIt != getEdgeIdx.end())
            {
                size_t edgeIdx = edgeIdxIt->second;
                currentStructuralAlignment.push_back(edgeIdx);
                inSolution[edgeIdx] = true;
                ++primal[edgeIdx]; // primal counter
            }
            else
            {
                std::cerr << "Alignment match where no alignment edge is defined!!" << std::endl;
                exit(1);
            }
        }
        dualValue = optScore;

        // store the subgradient values
        if (dualValue < bestUpperBoundScore)
        {
            for (size_t edgeIdx = 0ul; edgeIdx < sourceNode.size(); ++edgeIdx)
            {
                bestUpperBoundScores[sourceNode[edgeIdx]][targetNode[edgeIdx]] = maxProfit[edgeIdx];
                bestLagrangianMultipliers = dual;
                bestUpperBoundScore = dualValue;
            }
        }

        // store all upper bounds
        allUpperBounds.push_back(dualValue);
        subgradientIndices.clear();
        for (size_t idx : currentStructuralAlignment)
        {
            size_t const & maxPE = maxProfitEdge[idx];

            if (inSolution[maxPE] && maxProfitEdge[maxPE] == idx)
                continue;

            auto dualIt = pairedEdgesToDual.find(std::make_pair(idx, maxPE));
            if (dualIt != pairedEdgesToDual.end())
            {
                subgradient[dualIt->second] = 1.0;
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
                subgradient[dualIt->second] = -1.0;
                subgradientIndices.push_back(dualIt->second);
            }
            else
            {
                std::cerr << "Strange condition: pairedEdgesToDual(" << maxPE << "," << idx
                          << ") not defined!" << std::endl;
                exit(1);
            }
        }

        double lowerBound = 0.0;
        std::map<size_t, size_t> contacts;
        if (doMatching)
        {
            _VV(params, "start matching");
#ifdef LEMON_FOUND
            lowerBound = computeMatching(contacts, currentStructuralAlignment, inSolution);
#else
            std::cerr << "Cannot compute a matching without the Lemon library. Please install Lemon and try again."
                      << std::endl;
            exit(1);
#endif
        }
        else
        {
            _VV(params, "no matching");
            double seqPart = 0.0;
            double structPart = 0.0;

            std::map<PosPair, bool> computed;
            for (size_t line : currentStructuralAlignment)
            {
                seqPart += sequencesScore[line];
                size_t maxPE = maxProfitEdge[line];
                PosPair interaction{line, maxPE};
                PosPair revInteraction{maxPE, line};
                if (inSolution[maxPE] && line != maxPE && computed.count(interaction) == 0)
                {
                    if (structureScore[interaction] + structureScore[revInteraction] > 0)
                    {
                        structPart += structureScore[interaction] + structureScore[revInteraction];
                        computed[interaction] = true;
                        computed[revInteraction] = true;
                        _VV(params, "solution " << line << " " << maxPE << " " << structureScore[interaction] +
                                                                                 structureScore[revInteraction]);
                    }
                }
            }
            lowerBound = seqPart + structPart;
            _VV(params, "without matching: " << seqPart << " + " << structPart << " = " << lowerBound);
        }

        for (size_t idx : currentStructuralAlignment)
        {
            size_t const & maxPE = maxProfitEdge[idx];
            _VVV(params, "Alignment[" << idx << "; " << sourceNode[idx] << "," << targetNode[idx] << "] maxProfitEdge ["
                                      << maxPE << "; " << sourceNode[maxPE] << "," << targetNode[maxPE] << "] score "
                                      << maxProfit[idx] << " inSolution " << inSolution[maxPE] << " rec "
                                      << (maxProfitEdge[maxPE] == idx) << " mwm " << (contacts[idx] == maxPE));
        }

        // we have to substract the gapcosts, otherwise the lower bound might be higher than the upper bound
        primalValue = lowerBound + gapScore;
        _VV(params, "primal " << primalValue << " = " << lowerBound << " (lb) + " << gapScore << " (gap)");
        _VV(params, "dual " << scoreseq << " == " << dualValue);

        // store the best alignment found so far
        if (primalValue > bestStructuralAlignmentScore)
        {
            bestStructuralAlignmentScore = primalValue;
            bestStructuralAlignment = currentStructuralAlignment;
            edgeMatching = contacts;
            bestAlignment = Alignment(alignment);
        }
        ++numIterations;
    };

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
            structureLines.emplace_back(sourceNode[idx] + 1, targetNode[idx] + 1, (idx != edgeMatching.at(idx)));
        }
        return structureLines;
    }
};

} // namespace lara

