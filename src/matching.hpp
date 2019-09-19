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

/*!\file matching.hpp
 * \brief This file contains the maximum weighted matching algorithms.
 */

#include <cstdint>
#include <iostream>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef LEMON_FOUND
#include <lemon/core.h>
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/concepts/graph.h>
#endif

namespace lara
{

typedef std::tuple<float, size_t, size_t> Interaction;    // probability, lineL, lineR
typedef std::set<Interaction>::iterator InteractionIterator;
typedef std::pair<InteractionIterator, InteractionIterator> EdgeConflict;
typedef std::set<InteractionIterator> InteractionSet;

} // namespace lara

namespace std
{

template <>
struct less<lara::InteractionIterator>
{
    bool operator()(const lara::InteractionIterator & lhs, const lara::InteractionIterator & rhs) const
    {
        return make_pair(get<1>(*lhs), get<2>(*lhs)) < make_pair(get<1>(*rhs), get<2>(*rhs));
    }
};

} // namespace std

namespace lara
{

class Matching
{
private:
    std::unordered_map<size_t, size_t> contacts;
    std::vector<std::vector<Contact>> const & possiblePartners;
    size_t const algorithm;

    //!\brief Helper function that calculates whether two interactions use the same vertex.
    bool hasConflict(Interaction const & a, Interaction const & b)
    {
        return std::get<1>(a) == std::get<1>(b)
            || std::get<1>(a) == std::get<2>(b)
            || std::get<2>(a) == std::get<1>(b)
            || std::get<2>(a) == std::get<2>(b);
    }

    InteractionSet solveConflicts(float & weight, std::vector<EdgeConflict> const & conflicts)
    {
        if (conflicts.empty())
            return InteractionSet();

        auto edge_cmp = [] (InteractionIterator const & a, InteractionIterator const & b) { return (*a) >= (*b); };
        InteractionIterator edgeS = std::min(conflicts.front().first, conflicts.front().second, edge_cmp);
        InteractionIterator edgeL = std::max(conflicts.front().first, conflicts.front().second, edge_cmp);
        float weightS = -std::get<0>(*edgeS);
        float weightL = -std::get<0>(*edgeL);

        if (conflicts.size() == 1)
        {
            weight += weightS;
            return InteractionSet{edgeS};
        }

        std::vector<EdgeConflict> remainingS;
        std::copy_if(conflicts.begin(), conflicts.end(), std::back_inserter(remainingS),
                     [&edgeS] (EdgeConflict const & c) {return edgeS != c.first && edgeS != c.second;});
        InteractionSet eliminateS = solveConflicts(weightS, remainingS);

        if (weightS > weightL) // prune if S is already smaller
        {
            std::vector<EdgeConflict> remainingL;
            std::copy_if(conflicts.begin(), conflicts.end(), std::back_inserter(remainingL),
                         [&edgeL] (EdgeConflict const & c) {return edgeL != c.first && edgeL != c.second;});

            if (remainingS != remainingL) // prune if L has the same subtree
            {
                InteractionSet eliminateL = solveConflicts(weightL, remainingL);
                if (weightS > weightL)
                {
                    eliminateL.insert(edgeL);
                    weight += weightL;
                    return eliminateL;
                }
            }
        }
        eliminateS.insert(edgeS);
        weight += weightS;
        return eliminateS;
    }

    float computeGreedyMatching(std::vector<size_t> const & currentAlignment,
                                size_t lookahead = 5ul)
    {
        // fill the priority queue with interaction edges
        float score = 0.0f;
        std::set<Interaction> queue;
        for (size_t idx = 0ul; idx < currentAlignment.size(); ++idx)
            for (Contact const & contact : possiblePartners[idx])
                queue.emplace(-2 * contact.first, currentAlignment[idx], contact.second);

        if (lookahead > queue.size())
            lookahead = queue.size();
        else if (lookahead == 0)
            lookahead = 5;

        // select the best edges from queue
        std::vector<InteractionIterator> selection{};
        selection.reserve(lookahead);
        contacts.clear();
        auto queueIt = queue.begin();
        while (queueIt != queue.end())
        {
            for (selection.clear(); selection.size() < lookahead && queueIt != queue.end(); ++queueIt)
                if (contacts.count(std::get<1>(*queueIt)) == 0 && contacts.count(std::get<2>(*queueIt)) == 0)
                    selection.push_back(queueIt);

            // search conflicts
            std::vector<EdgeConflict> conflicts{};
            for (auto itA = selection.begin(); itA != selection.end(); ++itA)
                for (auto itB = itA + 1; itB != selection.end(); ++itB)
                    if (hasConflict(**itA, **itB))
                        conflicts.emplace_back(*itA, *itB);

            // solve conflicts
            float weight = 0.0f;
            InteractionSet eliminate = solveConflicts(weight, conflicts);

            // save MWM contacts and count score
            for (InteractionIterator const & it : selection)
            {
                if (eliminate.count(it) == 0)
                {
                    contacts[std::get<1>(*it)] = std::get<2>(*it);
                    contacts[std::get<2>(*it)] = std::get<1>(*it);
                    score += -std::get<0>(*it);
                }
            }
        }
        return score;
    }

#ifdef LEMON_FOUND
    /*!
     * \brief Computes a maximum weighted matching using LEMON.
     * \param[out] contacts         The pairing result of the matching.
     * \param[in]  currentAlignment The active lines, which build the alignment.
     * \returns The score of the matching.
     */
    float computeLemonMatching(std::vector<size_t> const & currentAlignment)
    {
        contacts.clear();

        lemon::SmartGraph lemonG;
        std::unordered_map<size_t, lemon::SmartGraph::Node> nodes{};

        for (size_t const & line : currentAlignment)
        {
            contacts[line] = line;
            nodes[line] = lemonG.addNode();
        }

        typedef lemon::SmartGraph::EdgeMap<float> EdgeMap;
        EdgeMap weight(lemonG);
        lemon::SmartGraph::EdgeMap<PosPair> interactions(lemonG);
        for (size_t idx = 0ul; idx < currentAlignment.size(); ++idx)
        {
            for (Contact const & contact : possiblePartners[idx])
            {
                auto newEdge = lemonG.addEdge(nodes[currentAlignment[idx]], nodes[contact.second]);
                weight[newEdge] = 2 * contact.first;
                interactions[newEdge] = std::make_pair(currentAlignment[idx], contact.second);
            }
        }
        lemon::MaxWeightedMatching<lemon::SmartGraph, EdgeMap> mwm(lemonG, weight);
        mwm.run();
        for (lemon::SmartGraph::EdgeIt edgeIt(lemonG); edgeIt!=lemon::INVALID; ++edgeIt)
        {
            PosPair inter = interactions[edgeIt];
            if (mwm.matching(edgeIt))
            {
                contacts[inter.first] = inter.second;
                contacts[inter.second] = inter.first;
            }
        }
        return mwm.matchingWeight();
    }
#endif

public:
    explicit Matching(std::vector<std::vector<Contact>> const & possiblePartners_, size_t algorithm_) :
        contacts(),
        possiblePartners(possiblePartners_),
        algorithm(algorithm_)
    {}

    std::unordered_map<size_t, size_t> getContacts()
    {
        return contacts;
    }

    float computeScore(std::vector<size_t> const & currentAlignment)
    {
#ifdef LEMON_FOUND
        if (algorithm == 0)
            return computeLemonMatching(currentAlignment);
        else
#endif
            return computeGreedyMatching(currentAlignment, algorithm);
    }
};

} // namespace lara
