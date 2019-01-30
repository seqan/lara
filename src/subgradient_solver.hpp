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

/*!\file subgradient_solver.hpp
 * \brief This file contains the Subgradient Solver for LaRA.
 */

#include <iostream>
#include <list>
#include <ostream>
#include <utility>
#include <vector>

#include "data_types.hpp"
#include "parameters.hpp"

namespace lara
{

class SubgradientSolver
{
public:
    Lagrange lagrange;
    float stepSizeFactor;
    float bestLowerBound;
    float bestUpperBound;
    float currentLowerBound;
    float currentUpperBound;
    size_t nondecreasingRounds;
    unsigned remainingIterations;
    PosPair sequenceIndices;

    std::vector<float> subgradient;
    std::vector<float> dual;
    std::list<size_t> subgradientIndices;

    SubgradientSolver(PosPair indices, InputStorage const & store, Parameters & params)
        : lagrange(store[indices.first], store[indices.second], params)
    {
        stepSizeFactor = params.stepSizeFactor;
        bestLowerBound = negInfinity;
        bestUpperBound = posInfinity;

        currentLowerBound = negInfinity;
        currentUpperBound = posInfinity;
        nondecreasingRounds = 0ul;
        sequenceIndices = indices;
        remainingIterations = params.numIterations;
        subgradient.resize(lagrange.getDimension().second);
        dual.resize(subgradient.size());
    }
};

class SubgradientSolverMulti
{
private:
    // parameters
    InputStorage const & store;
    Parameters & params;

    std::vector<PosPair> inputPairs;
    std::vector<SubgradientSolver> solvers;

public:
    SubgradientSolverMulti(InputStorage const & _store, Parameters & _params)
        : store(_store), params(_params)
    {
        inputPairs.reserve(store.size() * (store.size() - 1ul) / 2ul);
        for (size_t idxA = 0ul; idxA < store.size() - 1ul; ++idxA)
            for (size_t idxB = idxA + 1ul; idxB < store.size(); ++idxB)
                inputPairs.emplace_back(idxA, idxB);

        solvers.reserve(std::min((size_t)params.num_threads, inputPairs.size()));
    }

    float calcStepsize(SubgradientSolver const & slv) const
    {
        return slv.stepSizeFactor * (slv.bestUpperBound - slv.bestLowerBound) / slv.subgradientIndices.size();
    }

    float getLowerBound(uint8_t sIdx)
    {
        return solvers[sIdx].bestLowerBound;
    }

    float getUpperBound(uint8_t sIdx)
    {
        return solvers[sIdx].bestUpperBound;
    }

    void solve(lara::OutputTCoffeeLibrary & results)
    {
        std::vector<PosPair>::const_iterator inputPairIter;
        for (inputPairIter = inputPairs.cbegin();
             inputPairIter != inputPairs.cend() && solvers.size() < params.num_threads;
             ++inputPairIter)
        {
            solvers.emplace_back(*inputPairIter, store, params);
        }

        size_t num_at_work = solvers.size();
        std::vector<bool> at_work;
        at_work.resize(num_at_work, true);

        while (num_at_work > 0ul)
        {
            #pragma omp parallel for num_threads(params.num_threads)
            for (size_t idx = 0ul; idx < solvers.size(); ++idx)
            {
                if (!at_work[idx])
                    continue;

                solvers[idx].currentUpperBound = solvers[idx].lagrange.relaxed_solution(params.laraGapOpen,
                                                                                        params.laraGapExtend);

                solvers[idx].currentLowerBound = solvers[idx].lagrange.valid_solution(solvers[idx].subgradient,
                                                                                      solvers[idx].subgradientIndices);

//                _LOG(2, "(" << solvers[idx].remainingIterations << ") \tbest: " << solvers[idx].bestUpperBound << "\t/"
//                            << solvers[idx].bestLowerBound << "\t"
//                            << "current: " << solvers[idx].currentUpperBound << "/\t" << solvers[idx].currentLowerBound
//                            << "\t(" << solvers[idx].subgradientIndices.size() << ")\t");

                // compare upper and lower bound
                if (solvers[idx].currentUpperBound < solvers[idx].bestUpperBound)
                {
                    solvers[idx].bestUpperBound = solvers[idx].currentUpperBound;
                    solvers[idx].nondecreasingRounds = 0ul;
                }

                if (solvers[idx].currentLowerBound > solvers[idx].bestLowerBound)
                {
                    solvers[idx].bestLowerBound = solvers[idx].currentLowerBound;
                    solvers[idx].nondecreasingRounds = 0ul;
//                    _LOG(1, "(*)");

                }
//                _LOG(2, std::endl);

                if (solvers[idx].nondecreasingRounds++ >= params.maxNondecrIterations)
                {
                    solvers[idx].stepSizeFactor /= 2.0f;
                    solvers[idx].nondecreasingRounds = 0;
                }

                float stepSize = calcStepsize(solvers[idx]);
                for (size_t si : solvers[idx].subgradientIndices)
                {
                    solvers[idx].dual[si] -= stepSize * solvers[idx].subgradient[si];
                    solvers[idx].subgradient[si] = 0.0f;
                }
                --solvers[idx].remainingIterations;

                SEQAN_ASSERT_MSG(!solvers[idx].subgradientIndices.empty() ||
                                 solvers[idx].currentUpperBound - solvers[idx].currentLowerBound < 0.1f,
                                 "The bounds differ, although there are no subgradients.");
                SEQAN_ASSERT_GT_MSG(solvers[idx].bestUpperBound + params.epsilon, solvers[idx].bestLowerBound,
                                    "The lower boundary exceeds the upper boundary.");

                if (solvers[idx].bestUpperBound - solvers[idx].bestLowerBound < params.epsilon ||
                    solvers[idx].remainingIterations == 0u)
                {
                    #pragma omp critical
                    {
                        _LOG(1, "Thread " << idx << " finished alignment (" << solvers[idx].sequenceIndices.first
                                          << "," << solvers[idx].sequenceIndices.second << ")." << std::endl);
//                        lara::printAlignment(std::cerr,
//                                             solvers[idx].lagrange.getAlignment(),
//                                             store[solvers[idx].sequenceIndices.first].name,
//                                             store[solvers[idx].sequenceIndices.second].name);
                        results.addAlignment(solvers[idx].lagrange, solvers[idx].sequenceIndices);

                        if (inputPairIter == inputPairs.cend())
                        {
                            at_work[idx] = false;
                            --num_at_work;
                        }
                        else
                        {
                            solvers[idx] = SubgradientSolver(*inputPairIter, store, params);
                            ++inputPairIter;
                        }
                    }; // end critical region
                }
                else
                {
                    solvers[idx].lagrange.updateScores(solvers[idx].dual, solvers[idx].subgradientIndices);
                }
            }
        } // end while
    }
};

} // namespace lara

