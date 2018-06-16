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
    double epsilon;
    double stepSizeFactor;
    unsigned verbose;
    size_t numIterations;
    size_t maxNondecreasingIterations;
    size_t iterationIdx;

    std::vector<double> primal;
    std::vector<double> dual;
    std::vector<double> multiplierLowerBound;
    std::vector<double> multiplierUpperBound;

    double currentLowerBound;
    double currentUpperBound;
    double bestLowerBound;
    double bestUpperBound;

public:
    void initialize(PosPair const dim, Parameters const & params)
    {
        primal.resize(dim.first);
        dual.resize(dim.second);
        multiplierLowerBound.resize(dim.second, negInfinity);
        multiplierUpperBound.resize(dim.second, posInfinity);

        currentLowerBound = negInfinity;
        currentUpperBound = posInfinity;
        bestLowerBound    = negInfinity;
        bestUpperBound    = posInfinity;

        //        stepSize = 0.0;
        epsilon        = params.epsilon;
        stepSizeFactor = params.stepSizeFactor;
        verbose        = params.verbose;
        numIterations  = params.numIterations;
        maxNondecreasingIterations = params.numNondecreasingIterations;
        iterationIdx   = 0ul;
    }

    double calcStepsize(double upper, double lower, size_t numSubgradients)
    {
        return stepSizeFactor * (upper - lower) / numSubgradients;
    }

    double getLowerBound()
    {
        return bestLowerBound;
    }

    double getUpperBound()
    {
        return bestUpperBound;
    }

    Status solve(Lagrange & lagrange)
    {
        size_t nondecreasingRounds = 0ul;

        std::list<size_t> dualIndices;
        std::list<size_t> subgradientIndices;

        std::vector<double> subgradient;
        subgradient.resize(dual.size());

        for (iterationIdx = 0ul; iterationIdx < numIterations; ++iterationIdx)
        {
            lagrange.evaluate(dual, dualIndices, currentUpperBound, currentLowerBound, subgradient, subgradientIndices);

            if (verbose >= 1)
                std::cerr << "(" << iterationIdx << ") \tbest: " << bestUpperBound << "\t/" << bestLowerBound << "\t"
                          << "current: " << currentUpperBound << "/\t" << currentLowerBound
                          << "\t(" << subgradientIndices.size() << ")\t";

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
                if (verbose >= 1)
                    std::cerr << "(*)";

            }
            if (verbose >= 1)
                std::cerr << std::endl;

            if (bestUpperBound - bestLowerBound < epsilon)
                return Status::EXIT_OK;

            if (subgradientIndices.empty() && currentUpperBound - currentLowerBound > epsilon)
            {
                std::cerr << "Error: The bounds differ, although there are no subgradients." << std::endl;
                return Status::EXIT_ERROR;
            }

            if (nondecreasingRounds++ >= maxNondecreasingIterations)
            {
                stepSizeFactor /= 2.0;
                nondecreasingRounds = 0;
            }

            double stepSize = calcStepsize(bestUpperBound, bestLowerBound, subgradientIndices.size());
            if (verbose >= 2)
                std::cerr << "stepsize = " << stepSize << std::endl;

            bool changed = false;
            for (size_t idx : subgradientIndices)
            {
                double oldValue = dual[idx];
                double newValue = dual[idx] - stepSize * subgradient[idx];
                if (newValue < multiplierLowerBound[idx])
                    dual[idx] = multiplierLowerBound[idx];
                else if (newValue > multiplierUpperBound[idx])
                    dual[idx] = multiplierUpperBound[idx];
                else
                    dual[idx] = newValue;

                if (dual[idx] != oldValue)
                    changed = true;
            }
            if (!changed)
            {
                std::cerr << "Error: StepSizeFactor below precision." << std::endl;
                return Status::EXIT_ERROR;
            }

            if (verbose >= 2)
            {
                std::cerr << "iteration " << iterationIdx << std::endl;
                std::cerr << "subgradient.size  = " << subgradient.size() << std::endl;
                std::cerr << "multiplierLowerBound.size = " << multiplierLowerBound.size() << std::endl;
                std::cerr << "multiplierUpperBound.size = " << multiplierUpperBound.size() << std::endl;
            }

            dualIndices = std::list<size_t>(subgradientIndices);
            subgradientIndices.clear();

            for (size_t dualIdx : dualIndices)
            {
                subgradient[dualIdx] = 0.0;
            }
        } // end for (iterationIdx = 0..numIterations-1)

        return iterationIdx < numIterations ? Status::EXIT_OK : Status::CONTINUE;
    }
};

} // namespace lara

