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
//    double stepSize;
    double stepSizeFactor;
    size_t numIterations;
//    size_t maxNondecreasingIterations;
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
    SubgradientSolver(std::pair<size_t, size_t> const dim, Parameters const & params)
    {
//        stepSize = 0.0;
        stepSizeFactor = params.stepSizeFactor;
        numIterations = params.numIterations;
//        maxNondecreasingIterations = params.numNondecreasingIterations;
        iterationIdx = 0ul;

        primal.resize(dim.first);
        dual.resize(dim.second);
        multiplierLowerBound.resize(dim.second, negInfinity);
        multiplierUpperBound.resize(dim.second, posInfinity);

        currentLowerBound = negInfinity;
        currentUpperBound = posInfinity;
        bestLowerBound = negInfinity;
        bestUpperBound = posInfinity;
    }

    double calcStepsize()
    {
        return 1.0 / (iterationIdx + 1.0);
    }

    double calcStepsize(double upper, double lower, size_t numSubgradients)
    {
        return stepSizeFactor * (upper - lower) / numSubgradients;
    }

    Status solve(bool verbose, double epsilon, size_t maxNondecreasingIterations)
    {
        size_t nondecreasingRounds = 0ul;

        std::list<size_t> dualIndices;
        std::list<size_t> subgradientIndices;

        std::vector<double> subgradient;
        subgradient.resize(dual.size());

        std::vector<double> primalVector;
        primalVector.resize(primal.size());

        std::vector<double> primalFeasibleSolution;
        primalFeasibleSolution.resize(primal.size());

        for (iterationIdx = 0ul; iterationIdx < numIterations; ++iterationIdx)
        {
            if (verbose)
                std::cerr << "(" << iterationIdx << ")\tbest: " << bestUpperBound << "\t/" << bestLowerBound << "\t";

            // TODO lagrangeEvaluateProblem

            if (verbose)
                std::cerr << "current: " << currentUpperBound << "/\t" << currentLowerBound
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
                if (verbose)
                    std::cerr << "(*)";

            }
            if (verbose)
                std::cerr << std::endl;

            if (bestUpperBound - bestLowerBound < epsilon)
                return Status::EXIT_OK;

            if (subgradientIndices.empty() && currentUpperBound - currentLowerBound > epsilon)
            {
                std::cout << "Error: The bounds differ, although there are no subgradients." << std::endl;
                return Status::EXIT_ERROR;
            }

            if (nondecreasingRounds++ >= maxNondecreasingIterations)
            {
                stepSizeFactor /= 2.0;
                nondecreasingRounds = 0;
            }

            double stepSize = calcStepsize(bestUpperBound, bestLowerBound, subgradientIndices.size());
            if (verbose)
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
                std::cout << "Error: StepSizeFactor below precision." << std::endl;
                return Status::EXIT_ERROR;
            }

            if (verbose)
            {
                std::cout << "subgradient.size  = " << subgradient.size() << std::endl;
                std::cout << "primalVector.size = " << primalVector.size() << std::endl;
                std::cout << "primalfeasibleSolution.size = " << primalFeasibleSolution.size() << std::endl;
                std::cout << "multiplierLowerBound.size = " << multiplierLowerBound.size() << std::endl;
                std::cout << "multiplierUpperBound.size = " << multiplierUpperBound.size() << std::endl;
            }

            dualIndices = subgradientIndices;
            subgradientIndices.clear();

            for (size_t dualIdx : dualIndices)
            {
                subgradient[dualIdx] = 0.0;
            }
        } // end for (iterationIdx = 0..numIterations-1)

        return iterationIdx < numIterations ? Status::EXIT_OK : Status::CONTINUE;
    }
};

std::ostream & operator<<(std::ostream & stream, SubgradientSolver const & /* solver */)
{
    return stream;
}

} // namespace lara

