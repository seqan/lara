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
private:
    float epsilon;
    float stepSizeFactor;
    size_t numIterations;
    size_t maxNondecreasingIterations;

    std::vector<float> dual{};
    float currentLowerBound;
    float currentUpperBound;
    float bestLowerBound;
    float bestUpperBound;

public:
    SubgradientSolver(float eps, float stepFac, size_t nIter, size_t maxNondecrIter)
        : epsilon(eps), stepSizeFactor(stepFac), numIterations(nIter), maxNondecreasingIterations(maxNondecrIter)
    {
        currentLowerBound = negInfinity;
        currentUpperBound = posInfinity;
        bestLowerBound    = negInfinity;
        bestUpperBound    = posInfinity;
    }

    float calcStepsize(float upper, float lower, size_t numSubgradients)
    {
        return stepSizeFactor * (upper - lower) / numSubgradients;
    }

    float getLowerBound()
    {
        return bestLowerBound;
    }

    float getUpperBound()
    {
        return bestUpperBound;
    }

    Status solve(Lagrange & lagrange)
    {
        size_t iterationIdx;
        size_t nondecreasingRounds = 0ul;

        std::list<size_t> subgradientIndices;
        std::vector<float> subgradient;

        dual.resize(lagrange.getDimension().second);
        subgradient.resize(dual.size());

        for (iterationIdx = 0ul; iterationIdx < numIterations; ++iterationIdx)
        {
            lagrange.updateScores(dual, subgradientIndices);
            currentUpperBound = lagrange.relaxed_solution();
            currentLowerBound = lagrange.valid_solution(subgradient, subgradientIndices);

            _LOG(1, "(" << iterationIdx << ") \tbest: " << bestUpperBound << "\t/" << bestLowerBound << "\t"
                        << "current: " << currentUpperBound << "/\t" << currentLowerBound
                        << "\t(" << subgradientIndices.size() << ")\t");

            // compare upper and lower bound
            if (currentUpperBound < bestUpperBound)
            {
                bestUpperBound      = currentUpperBound;
                nondecreasingRounds = 0ul;
            }

            if (currentLowerBound > bestLowerBound)
            {
                bestLowerBound      = currentLowerBound;
                nondecreasingRounds = 0ul;
                _LOG(1, "(*)");

            }
            _LOG(1, std::endl);

            if (subgradientIndices.empty() && currentUpperBound - currentLowerBound > 0.1f)
            {
                std::cerr << "Error: The bounds differ, although there are no subgradients." << std::endl;
                return Status::EXIT_ERROR;
            }

            if (bestUpperBound - bestLowerBound < epsilon)
            {
                SEQAN_ASSERT_GT(bestUpperBound + epsilon, bestLowerBound);
                _LOG(1, "Found the optimal alignment in iteration " << iterationIdx << "." << std::endl);
                return Status::EXIT_OK;
            }

            if (nondecreasingRounds++ >= maxNondecreasingIterations)
            {
                stepSizeFactor /= 2.0f;
                nondecreasingRounds = 0;
            }

            float stepSize = calcStepsize(bestUpperBound, bestLowerBound, subgradientIndices.size());
            _LOG(2, "stepsize = " << stepSize << std::endl);

            for (size_t idx : subgradientIndices)
            {
                dual[idx] -= stepSize * subgradient[idx];
                subgradient[idx] = 0.0f;
            }
        } // end for (iterationIdx = 0..numIterations-1)

        return iterationIdx < numIterations ? Status::EXIT_OK : Status::CONTINUE;
    }
};

} // namespace lara

