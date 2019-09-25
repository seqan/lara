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

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include <chrono>
#include <iostream>
#include <iterator>
#include <list>
#include <ostream>
#include <utility>
#include <vector>

#include <seqan/sequence.h>
#include <seqan/simd.h>

#include "data_types.hpp"
#include "lagrange.hpp"
#include "parameters.hpp"
#include "score.hpp"

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

    std::vector<float> subgradient{};
    std::vector<float> dual{};
    std::list<size_t> subgradientIndices{};

    SubgradientSolver(PosPair indices,
                      InputStorage const & store,
                      Parameters & params,
                      RnaScoreType * score,
                      size_t seqIdx):
        lagrange(store[indices.first], store[indices.second], params, score, seqIdx),
        stepSizeFactor{params.stepSizeFactor},
        bestLowerBound{negInfinity},
        bestUpperBound{posInfinity},
        currentLowerBound{negInfinity},
        currentUpperBound{posInfinity},
        nondecreasingRounds{0ul},
        remainingIterations{params.numIterations},
        sequenceIndices{indices}
    {
        subgradient.resize(lagrange.getDimension());
        dual.resize(subgradient.size());
    }

    SubgradientSolver()                                      = delete;
    SubgradientSolver(SubgradientSolver const &)             = default;
    SubgradientSolver(SubgradientSolver &&)                  = default;
    SubgradientSolver & operator=(SubgradientSolver const &) = default;
    SubgradientSolver & operator=(SubgradientSolver &&)      = default;
    ~SubgradientSolver()                                     = default;
};

class SubgradientSolverMulti
{
private:
    // parameters
    InputStorage const & store;
    Parameters & params;

    // sort alignment order for the first sequence length
    struct longer_seq
    {
        InputStorage const & store;

        bool operator()(PosPair lhs, PosPair rhs) const
        {
            if (seqan::length(store[lhs.first].sequence) > seqan::length(store[rhs.first].sequence))
                return true;
            else if (seqan::length(store[lhs.first].sequence) == seqan::length(store[rhs.first].sequence))
                return seqan::length(store[lhs.second].sequence) >= seqan::length(store[rhs.second].sequence);
            else
                return false;
        }
    };

    std::set<PosPair, longer_seq> inputPairs;
    std::vector<SubgradientSolver> solvers;

public:
    explicit SubgradientSolverMulti(InputStorage const & _store, Parameters & _params)
        : store(_store), params(_params), inputPairs(longer_seq{store})
    {
        // Add the sequence index pairs to the set of alignments, longer sequence first.
        for (size_t idxA = 0ul; idxA < store.size() - 1ul; ++idxA)
            for (size_t idxB = idxA + 1ul; idxB < store.size(); ++idxB)
                if (seqan::length(store[idxA].sequence) >= seqan::length(store[idxB].sequence))
                    inputPairs.emplace(idxA, idxB);
                else
                    inputPairs.emplace(idxB, idxA);
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
        _LOG(1, "3) Solve " << inputPairs.size() << " structural alignments..." << std::endl);
        if (inputPairs.empty())
            return;

#ifdef SEQAN_SIMD_ENABLED
        size_t const simd_len = seqan::LENGTH<typename seqan::SimdVector<ScoreType>::Type>::VALUE;
#else
        size_t const simd_len = 1ul;
#endif

        // Determine number of parallel alignments.
        Clock::time_point timeInit = Clock::now();
        size_t const num_parallel = std::min(simd_len * params.threads, inputPairs.size());
        size_t const num_threads = (num_parallel - 1) / simd_len + 1;

        // We iterate over all pairs of input sequences, starting with the longest.
        auto iter = inputPairs.cbegin(); // iter -> pair of sequence indices

        // Initialise the alignments.
        std::vector<std::pair<seqan::StringSet<GappedSeq>, seqan::StringSet<GappedSeq>>> alignments(num_threads);

        // Integer sequence from 0 until length of longest seq -1
        seqan::String<unsigned> integerSeq;
        seqan::resize(integerSeq, seqan::length(store[iter->first].sequence));
        std::iota(begin(integerSeq), end(integerSeq), 0u);

        // Store the integer sequences.
        using PrefixType = seqan::Prefix<seqan::String<unsigned>>::Type;
        seqan::StringSet<seqan::String<unsigned>> seq1;
        seqan::StringSet<seqan::String<unsigned>> seq2;
        seqan::reserve(seq1, num_parallel);
        seqan::reserve(seq2, num_parallel);

        // Initialise the scores.
        std::vector<RnaScoreType> scores(num_threads);
        size_t const max_2nd_length = seqan::length(store[iter->second].sequence);
        auto const go = static_cast<ScoreType>(params.rnaScore.data_gap_open * factor2int);
        auto const ge = static_cast<ScoreType>(params.rnaScore.data_gap_extend * factor2int);

        // Initialise the solvers.
        solvers.reserve(num_parallel);

        for (iter = inputPairs.cbegin(); solvers.size() < num_parallel; ++iter)
        {
            size_t const aliIdx = solvers.size() / simd_len;
            size_t const seqIdx = solvers.size() % simd_len;
            auto const len = std::make_pair(length(store[iter->first].sequence), length(store[iter->second].sequence));

            // Once for each chunk of size simd_len.
            if (seqIdx == 0)
            {
                // Initialise the alignments.
                seqan::reserve(alignments[aliIdx].first, simd_len);
                seqan::reserve(alignments[aliIdx].second, simd_len);

                // Initialise the scores.
                scores[aliIdx].init(len.first, std::min(max_2nd_length, len.first), go, ge);
                _LOG(2, "     Resize matrix: " << len.first << "*" << std::min(max_2nd_length, len.first) << std::endl);
            }

            // Fill the alignments.
            appendValue(seq1, PrefixType(integerSeq, len.first));
            appendValue(seq2, PrefixType(integerSeq, len.second));
            appendValue(alignments[aliIdx].first, GappedSeq(seqan::back(seq1)));
            appendValue(alignments[aliIdx].second, GappedSeq(seqan::back(seq2)));

            // Fill the solvers.
            solvers.emplace_back(*iter, store, params, &(scores[aliIdx]), seqIdx);
        }
        SEQAN_ASSERT_EQ(num_parallel, solvers.size());
        SEQAN_ASSERT_EQ(num_parallel, seqan::length(seq1));
        SEQAN_ASSERT_EQ(num_parallel, seqan::length(seq2));
        SEQAN_ASSERT_EQ(num_threads, seqan::length(alignments));
        _LOG(1, "   * set up initial " << num_parallel << " structural alignments -> " << timeDiff(timeInit) << "ms"
                << std::endl);

        Clock::duration durationAlign{};
        Clock::duration durationMatching{};
        Clock::duration durationUpdate{};
        Clock::duration durationSerial{};
        Clock::time_point timeIter = Clock::now();

        // in parallel for each (SIMD) alignment
        #pragma omp parallel for num_threads(params.threads)
        for (size_t aliIdx = 0ul; aliIdx < num_threads; ++aliIdx)
        {
            Clock::duration durationThreadAlign{};
            Clock::duration durationThreadMatching{};
            Clock::duration durationThreadUpdate{};
            Clock::time_point timeThreadSerial = Clock::now();
            size_t num_at_work = seqan::length(alignments[aliIdx].first);
            std::vector<bool> at_work(num_at_work, true);
            auto const interval = std::make_pair(aliIdx * simd_len, std::min((aliIdx + 1) * simd_len, num_parallel));
            scores[aliIdx].updateLongestSeq(seq1, seq2, interval);

            // loop the thread until there is no more work to do
            while (num_at_work > 0ul)
            {
                // Performs the structural alignment. Returns the dual value (upper bound, solution of relaxed problem).
                Clock::time_point timeCurrent = Clock::now();
                seqan::String<ScoreType> res = seqan::globalAlignment(alignments[aliIdx].first,
                                                                      alignments[aliIdx].second,
                                                                      scores[aliIdx]);
                durationThreadAlign += Clock::now() - timeCurrent;

                // Evaluate each alignment result and adapt multipliers.
                for (size_t idx = interval.first; idx < interval.second; ++idx)
                {
                    size_t const seqIdx = idx % simd_len;
                    if (!at_work[seqIdx])
                        continue;

                    SubgradientSolver & ss = solvers[idx];
                    ss.currentUpperBound = res[seqIdx] / factor2int; // global alignment result

                    timeCurrent = Clock::now();
                    ss.currentLowerBound = ss.lagrange.valid_solution(ss.subgradient, ss.subgradientIndices,
                                                                      std::make_pair(alignments[aliIdx].first[seqIdx],
                                                                                     alignments[aliIdx].second[seqIdx]),
                                                                      params.matching,
                                                                      params.rnaScore);
                    durationThreadMatching += Clock::now() - timeCurrent;

                    // compare upper and lower bound
                    if (ss.currentUpperBound < ss.bestUpperBound)
                    {
                        ss.bestUpperBound = ss.currentUpperBound;
                        ss.nondecreasingRounds = 0ul;
                    }

                    if (ss.currentLowerBound > ss.bestLowerBound)
                    {
                        ss.bestLowerBound = ss.currentLowerBound;
                        ss.nondecreasingRounds = 0ul;

                    }

                    if (ss.nondecreasingRounds++ >= params.maxNondecrIterations)
                    {
                        ss.stepSizeFactor /= 2.0f;
                        ss.nondecreasingRounds = 0;
                    }

                    float stepSize = calcStepsize(ss);
                    for (size_t si : ss.subgradientIndices)
                    {
                        ss.dual[si] -= stepSize * ss.subgradient[si];
                        ss.subgradient[si] = 0.0f;
                    }
                    --ss.remainingIterations;

                    SEQAN_ASSERT_MSG(!ss.subgradientIndices.empty() ||
                                     ss.currentUpperBound - ss.currentLowerBound < 0.1f,
                                     (std::string{"The bounds differ, although there are no subgradients. "} +
                                         "Problem in aligning sequences " +
                                         seqan::toCString(store[ss.sequenceIndices.first].name) + " and " +
                                         seqan::toCString(store[ss.sequenceIndices.second].name)).c_str());
                    SEQAN_ASSERT_GT_MSG(ss.bestUpperBound + params.epsilon,
                                        ss.bestLowerBound,
                                        (std::string{"The lower boundary exceeds the upper boundary. "} +
                                            "Problem in aligning sequences " +
                                            seqan::toCString(store[ss.sequenceIndices.first].name) + " and " +
                                            seqan::toCString(store[ss.sequenceIndices.second].name)).c_str());

                    // The alignment is finished.
                    if (ss.bestUpperBound - ss.bestLowerBound < params.epsilon || ss.remainingIterations == 0u)
                    {
                        PosPair currentSeqIdx{};
                        #pragma omp critical (finished_alignment)
                        {
                            // write results
                            results.addAlignment(ss.lagrange, ss.sequenceIndices, params);
                            _LOG(2, "     Thread " << aliIdx << "." << seqIdx << " finished alignment "
                                    << ss.sequenceIndices.first << "/" << ss.sequenceIndices.second << std::endl);

                            if (iter == inputPairs.cend())
                            {
                                at_work[seqIdx] = false;
                                --num_at_work;
                            }
                            else
                            {
                                currentSeqIdx = *iter;
                                ++iter;
                            }
                        } // end critical region

                        if (at_work[seqIdx])
                        {
                            // Reset scores.
                            scores[aliIdx].reset(seqIdx);

                            // Set new sequences.
                            seq1[idx] = PrefixType(integerSeq, length(store[currentSeqIdx.first].sequence));
                            seq2[idx] = PrefixType(integerSeq, length(store[currentSeqIdx.second].sequence));
                            alignments[aliIdx].first[seqIdx] = GappedSeq(seq1[idx]);
                            alignments[aliIdx].second[seqIdx] = GappedSeq(seq2[idx]);

                            // Set new score matrix.
                            solvers[idx] = SubgradientSolver(currentSeqIdx, store, params, &(scores[aliIdx]), seqIdx);
                            scores[aliIdx].updateLongestSeq(seq1, seq2, interval);
                        }
                    }
                    else
                    {
                        timeCurrent = Clock::now();
                        ss.lagrange.updateScores(ss.dual, ss.subgradientIndices, params.rnaScore);
                        durationThreadUpdate += Clock::now() - timeCurrent;
                    }
                }
            }

            #pragma omp critical (update_time)
            {
                durationAlign += durationThreadAlign;
                durationMatching += durationThreadMatching;
                durationUpdate += durationThreadUpdate;
                durationSerial += Clock::now() - timeThreadSerial;
            }
        } // end parallel for

        auto durationToSeconds = [] (Clock::duration duration)
            { return std::chrono::duration_cast<std::chrono::seconds>(duration).count(); };

        _LOG(1, "   * parallel iterations -> " << timeDiff<std::chrono::seconds>(timeIter) << "s" << std::endl
                << "     (serial: " << durationToSeconds(durationSerial)
                << "s, ali: " << durationToSeconds(durationAlign)
                << "s, match: " << durationToSeconds(durationMatching)
                << "s, update: " << durationToSeconds(durationUpdate) << "s)" << std::endl);
    }
};

} // namespace lara
