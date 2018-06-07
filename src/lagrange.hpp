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

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include <seqan/rna_io.h>
#include <seqan/score.h>

#include "data_types.hpp"
#include "score.hpp"

namespace lara
{

class Lagrange
{
private:
    Parameters const & params;

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

    // scoring system holding the scores for different pairs of residues
//    RnaScoreMatrix score = params.laraScoreMatrix;

    // number of iterations performed so far
    size_t numIterations;

    // number of residues in total
    size_t residueCount;

    // best upper bound found so far
    double bestUpperBound;
    double lastUpperBound;
    std::vector<double> allUpperBounds;

    // best alignment found so far as aligned strings
    std::string alignedSeqA;
    std::string alignedSeqB;

    seqan::Rna5String sequenceA;
    seqan::Rna5String sequenceB;


    bool stacking;
    bool doMatching; // True for bpp matrix input

    // the best structural alignment score found so far
    double bestStructuralAlignmentScore;

    // the best (lowest) upper bound found so far
    double bestUpperBoundScore;

    // number of the conserved base pairs in the best solution
    size_t bestAlignNumConservedBP;

    std::vector<std::pair<seqan::Rna5, seqan::Rna5>> alignmentEdges;
    std::vector<size_t> activeEdges;
    std::vector<size_t> sourceNode;
    std::vector<size_t> targetNode;
    std::map<PosPair, size_t> nodeToEdge;

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
//            _VV(params, "base pair " << origin << "\t" << partner << "\t" << probability);
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
    double evaluateLines(std::vector<PosPair> & lines, AlignmentRow const & rowA, AlignmentRow const & rowB)
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

public:
    Lagrange(seqan::RnaRecord const & recordA, seqan::RnaRecord const & recordB, Parameters const & _params)
        : params(_params)
    {
        _VV(params, recordA.sequence << std::endl << recordB.sequence);
        PosPair seqLen{seqan::length(recordA.sequence), seqan::length(recordB.sequence)};

        sequenceA = seqan::Rna5String{recordA.sequence};
        sequenceB = seqan::Rna5String{recordB.sequence};

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
        alignedSeqA = "";
        alignedSeqB = "";
        dimension.first = 0ul;
        dimension.second = 0ul;
        stacking = false;
        doMatching = false;

        // score of best alignment
        bestStructuralAlignmentScore = negInfinity;
        bestUpperBoundScore = posInfinity;
        bestAlignNumConservedBP = 0ul;

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

        // skip edge filter for now
        size_t numEdges = seqLen.first * seqLen.second;
        std::vector<PosPair> edges;
        edges.reserve(numEdges);
        for (size_t idxA = 0ul; idxA < seqLen.first; ++idxA)
        {
            for (size_t idxB = 0ul; idxB < seqLen.second; ++idxB)
            {
                edges.emplace_back(idxA, idxB);
            }
        }
        _VV(params, "Created " << edges.size() << " = " << numEdges << " edges.");

        alignmentEdges.reserve(numEdges);
        activeEdges.reserve(numEdges);
        sourceNode.reserve(numEdges);
        targetNode.reserve(numEdges);
        primal.resize(numEdges, 0.0);
        sequencesScore.resize(numEdges, 0.0);
        possiblePartners.resize(numEdges);
        maxProfit.resize(numEdges, negInfinity);
        maxProfitEdge.resize(numEdges, SIZE_MAX);
        priorityQ.resize(numEdges);
        dualToPairedEdges.reserve(numEdges);

        size_t edgeIdx = 0ul;
        for (PosPair & edge : edges)
        {
            alignmentEdges.emplace_back(recordA.sequence[edge.first], recordB.sequence[edge.second]);
            activeEdges.push_back(edgeIdx);
            sourceNode.push_back(edge.first);
            targetNode.push_back(edge.second);
            nodeToEdge[edge] = edgeIdx++;
//            _VV(params, "edge " << edge.first << " " << edge.second);
        }
        bppGraphs = std::make_pair(seqan::front(recordA.bppMatrGraphs), seqan::front(recordB.bppMatrGraphs));
    }

    void start()
    {
        size_t dualIdx = 0ul;
        size_t numPartners = 0ul;
        for (size_t edgeIdx : activeEdges)
        {
            size_t headNode = sourceNode[edgeIdx];
            size_t tailNode = targetNode[edgeIdx];
            std::vector<Contact> headContact;
            std::vector<Contact> tailContact;
            extractContacts(headContact, bppGraphs.first, headNode);
            extractContacts(tailContact, bppGraphs.second, tailNode);
//            _VV(params, "active edge " << edgeIdx << " (" << headNode << "," << tailNode << ")");

            double alignScore = seqan::score(params.laraScoreMatrix,
                                             alignmentEdges[edgeIdx].first,
                                             alignmentEdges[edgeIdx].second);
//            _VV(params, "align score = " << alignScore);

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
                        auto partnerIter = nodeToEdge.find(std::make_pair(head.first, tail.first));
                        if (partnerIter != nodeToEdge.end() && nonCrossingEdges(edgeIdx, partnerIter->second))
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
//                            _VV(params, "dual idx " << dualIdx << " = (" << sourceNode[edgeIdx]
//                                                    << "-" << targetNode[edgeIdx] << ") -> (" << sourceNode[partnerIdx]
//                                                    << "-" << targetNode[partnerIdx] << ")");
                            pairedEdgesToDual[interaction] = dualIdx++;
                            dualToPairedEdges.emplace_back(edgeIdx, partnerIdx);

                            if (structScore + alignScore > maxProfit[edgeIdx])
                            {
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
        _VV(params, "Average number of partner edges = " << 1.0 * numPartners / activeEdges.size());
        dimension.second = dualIdx;
        dimension.first = activeEdges.size();
    }

    void evaluate(std::vector<double> & dual,
                  std::list<size_t> & dualIndices,
                  double & dualValue,   // upper bound
                  double & primalValue, // lower bound
                  std::vector<double> & subgradient,
                  std::list<size_t> & subgradientIndices,
                  std::vector<double> & primalSolution)
    {
        _VV(params, "evaluate()");
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
        for (size_t edgeIdx : activeEdges)
        {
            maxProfitScores[sourceNode[edgeIdx]][targetNode[edgeIdx]] = maxProfit[edgeIdx];
        }

        seqan::Score<double, seqan::RnaStructureScore> scoreAdaptor(&maxProfitScores,
                                                                    params.laraGapOpen,
                                                                    params.laraGapExtend);
        Alignment alignment;
        seqan::resize(seqan::rows(alignment), 2);
        AlignmentRow & rowA = seqan::row(alignment, 0);
        AlignmentRow & rowB = seqan::row(alignment, 1);
        std::vector<PosPair> lines{};
        _VV(params, "sequence length " << length(sequenceA) << " " << length(sequenceB));

        seqan::assignSource(rowA, sequenceA);
        seqan::assignSource(rowB, sequenceB);

        double optScore = seqan::globalAlignment(alignment, scoreAdaptor, seqan::AffineGaps());
        double gapScore = evaluateLines(lines, rowA, rowB);
        std::cerr << "Score: " << optScore << std::endl << alignment << std::endl;

        std::list<size_t> partnerEdgeIndices;
        std::map<size_t, bool> inSolutionMap;
        std::vector<bool> inSolutionVector;
        inSolutionVector.resize(alignmentEdges.size(), false);

        for (PosPair line : lines)
        {
            auto partnerIt = nodeToEdge.find(std::make_pair(line.first, line.second));
            if (partnerIt != nodeToEdge.end())
            {
                partnerEdgeIndices.push_back(partnerIt->second);
                inSolutionMap[partnerIt->second] = true;
                inSolutionVector[partnerIt->second] = true;
                ++primal[partnerIt->second]; // primal counter
                if (primalSolution.size() > 0)
                    primalSolution[partnerIt->second] = 1.0; // TODO check this
            }
            else
            {
                std::cerr << "Alignment match where no alignment edge is defined!!" << std::endl;
            }
        }
        dualValue = optScore;

        // store the subgradient values
        if (dualValue < bestUpperBoundScore)
        {
            for (size_t edgeIdx : activeEdges)
            {
                bestUpperBoundScores[sourceNode[edgeIdx]][targetNode[edgeIdx]] = maxProfit[edgeIdx];
                bestLagrangianMultipliers = dual;
                bestUpperBoundScore = dualValue;
            }
        }

        // store all upper bounds
        allUpperBounds.push_back(dualValue);
        lastUpperBound = dualValue;

        subgradientIndices.clear();
        for (size_t idx : partnerEdgeIndices)
        {
            size_t const maxPE = maxProfitEdge[idx];
            if (inSolutionVector[maxPE] && maxProfitEdge[maxPE] == idx)
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
            }
        }

        double lowerBound = 0.0;
        size_t numConserved = 0ul;
        if (doMatching)
        {
            _VV(params, "start matching");
            // TODO calculate lowerBound and numConserved
        }
        else
        {
            _VV(params, "no matching");
            // TODO calculate lowerBound
        }

        // we have to substract the gapcosts, otherwise the lower bound might be higher than the upper bound
        primalValue = lowerBound + gapScore;

        // store the best alignment found so far
        if (primalValue > bestStructuralAlignmentScore)
        {
            bestStructuralAlignmentScore = primalValue;
            bestAlignNumConservedBP = numConserved;
            // TODO bestStructural...
            // TODO alignedSeq
        }
        ++numIterations;
    };

    PosPair getDimension()
    {
        return dimension;
    }
};

//std::ostream & operator<<(std::ostream & stream, Lagrange const & lagrange)
//{
//    return stream << "Alignment" << std::endl;
//}

} // namespace lara

