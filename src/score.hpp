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

/*!\file alignment.hpp
 * \brief This file contains a position-dependent RNA scoring matrix.
 */

#include <iostream>
#include <iomanip>
#include <ostream>
#include <sstream>

#include <seqan/consensus.h>
#include <seqan/score.h>

namespace seqan
{

// --------------------------------------------------------
// POSITION DEPENDENT SCORE
// --------------------------------------------------------

struct RnaStructureScore;

template <typename ScoreType>
class Score<ScoreType, RnaStructureScore>
{
public:
    String<ScoreType, Alloc<OverAligned>> matrix; // aligned alloc
    ScoreType data_gap_open;
    ScoreType data_gap_extend;
    size_t dim; // length of second sequence
};

template <typename ScoreType, typename TPos>
inline
ScoreType score(Score<ScoreType, RnaStructureScore> const & sc, TPos const entryH, TPos const entryV)
{
    return sc.matrix[sc.dim * entryH + entryV];
}

#ifdef SEQAN_SIMD_ENABLED

// --------------------------------------------------------
// POSITION DEPENDENT SCORE WITH SIMD
// --------------------------------------------------------

struct RnaStructureScoreSimd;

template <typename ScoreType>
class Score<ScoreType, RnaStructureScoreSimd>
{
public:
    String<typename seqan::SimdVector<ScoreType>::Type, Alloc<OverAligned>> matrix; // aligned alloc
    ScoreType data_gap_open;
    ScoreType data_gap_extend;
    size_t dim; // length of second sequence
};

template <typename ScoreType, typename TPosVec>
inline
typename seqan::SimdVector<ScoreType>::Type score(Score<ScoreType, RnaStructureScoreSimd> const & sc,
                                                  TPosVec const entryH,
                                                  TPosVec const entryV)
{
    auto tmp = createVector<TPosVec>(sc.dim) * entryH + entryV;

    typename seqan::SimdVector<ScoreType>::Type res{};
    for (auto idx = 0u; idx < LENGTH<TPosVec>::VALUE; ++idx)
        res[idx] = sc.matrix[tmp[idx]][idx];

    return res;
}

// --------------------------------------------------------
// WRAPPER FOR POSITION DEPENDENT SCORE WITH SIMD
// --------------------------------------------------------

template <typename TScoreVec>
class Score<TScoreVec, ScoreSimdWrapper<Score<int32_t, RnaStructureScoreSimd>>>
{
public:
    using TVecValue = typename Value<TScoreVec>::Type;
    using TBaseScoreSpec = typename Spec<Score<int32_t, RnaStructureScoreSimd>>::Type;
    using TBaseScore = Score<typename IfC<sizeof(TVecValue) <= 2,
                                          int32_t,
                                          typename IfC<sizeof(TVecValue) == 8, int64_t, TVecValue>::Type
                                         >::Type,
                             TBaseScoreSpec>;

    // The score type is a ScoreMatrix.
    TScoreVec data_gap_extend = createVector<TScoreVec>(-1);
    TScoreVec data_gap_open   = createVector<TScoreVec>(-1);
    TBaseScore _baseScore;

    // Default Constructor.
    Score()
    {}

    template <typename TScoreVal2, typename TScoreSpec2>
    Score(Score<TScoreVal2, TScoreSpec2> const & pScore) :
        data_gap_extend(createVector<TScoreVec>(scoreGapExtend(pScore))),
        data_gap_open(createVector<TScoreVec>(scoreGapOpen(pScore))),
        _baseScore(pScore)
    {}
};

template <typename TValue, typename TVal1, typename TVal2>
TValue score(Score<TValue, ScoreSimdWrapper<Score<int32_t, RnaStructureScoreSimd>>> const & sc,
             TVal1 const & val1,
             TVal2 const & val2)
{
    return score(sc._baseScore, val1, val2);
}

#endif

} // namespace seqan
