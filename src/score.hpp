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

/*!\file score.hpp
 * \brief This file contains a position-dependent scoring matrix.
 */

#include <algorithm> // std::fill
#include <limits>    // std::numeric_limits

#include <seqan/score.h>
#include <seqan/sequence.h>
#include <seqan/simd.h>

namespace seqan
{

// --------------------------------------------------------
// POSITION SPECIFIC SCORE
// --------------------------------------------------------

struct PositionSpecificScore;

template <typename TScore>
class Score<TScore, PositionSpecificScore>
{
private:
    static TScore const INITVALUE;

public:
    String<TScore, Alloc<OverAligned>> matrix; // aligned alloc
    size_t dim; // length of second sequence
    TScore data_gap_open;
    TScore data_gap_extend;

    void init(size_t dim1, size_t dim2, TScore gapOpen, TScore gapExtend)
    {
        resize(matrix, dim1 * dim2, INITVALUE);
        dim = dim2;
        data_gap_open = gapOpen;
        data_gap_extend = gapExtend;
    }

    void set(size_t /* unused */, size_t idx1, size_t idx2, TScore value)
    {
        matrix[dim * idx1 + idx2] = value;
    }

    void reset(size_t /* unused */)
    {
        std::fill(begin(matrix), end(matrix), INITVALUE);
    }

    void updateLongestSeq(uint8_t /* unused */, uint8_t /* unused */) {}
};

template <typename TScore>
TScore const Score<TScore, PositionSpecificScore>::INITVALUE = std::numeric_limits<TScore>::lowest() / 3 * 2;

template <typename TScore, typename TPos>
inline
TScore score(Score<TScore, PositionSpecificScore> const & sc, TPos const entryH, TPos const entryV)
{
    return sc.matrix[sc.dim * entryH + entryV];
}

#ifdef SEQAN_SIMD_ENABLED

// --------------------------------------------------------
// POSITION SPECIFIC SCORE WITH SIMD
// --------------------------------------------------------

struct PositionSpecificScoreSimd;

template <typename TScore>
class Score<TScore, PositionSpecificScoreSimd>
{
private:
    using SimdScoreType = typename seqan::SimdVector<TScore>::Type;
    static TScore const INITVALUE;

public:
    String<SimdScoreType, Alloc<OverAligned>> matrix; // aligned alloc
    size_t dim; // length of second sequence
    TScore data_gap_open;
    TScore data_gap_extend;
    uint8_t longestSeqIdx1;
    uint8_t longestSeqIdx2;

    void init(size_t dim1, size_t dim2, TScore gapOpen, TScore gapExtend)
    {
        resize(matrix, dim1 * dim2, createVector<SimdScoreType>(INITVALUE));
        dim = dim2;
        data_gap_open = gapOpen;
        data_gap_extend = gapExtend;
        longestSeqIdx1 = 0;
        longestSeqIdx2 = 0;
    }

    void set(size_t seq, size_t idx1, size_t idx2, TScore value)
    {
        matrix[dim * idx1 + idx2][seq] = value;
    }

    void reset(size_t seq)
    {
        for (auto & score : matrix)
            score[seq] = INITVALUE;
    }

    void updateLongestSeq(uint8_t idx1, uint8_t idx2)
    {
        longestSeqIdx1 = idx1;
        longestSeqIdx2 = idx2;
    }
};

template <typename TScore>
TScore const Score<TScore, PositionSpecificScoreSimd>::INITVALUE = std::numeric_limits<TScore>::lowest() / 3 * 2;

template <typename TScore, typename TPosVec>
inline
typename seqan::SimdVector<TScore>::Type score(Score<TScore, PositionSpecificScoreSimd> const & sc,
                                                  TPosVec const entryH,
                                                  TPosVec const entryV)
{
    // The entry vectors consist of the same values, so we can always use the entry of the longest sequence.
    return sc.matrix[sc.dim * entryH[sc.longestSeqIdx1] + entryV[sc.longestSeqIdx2]];
}

// --------------------------------------------------------
// WRAPPER FOR POSITION SPECIFIC SCORE WITH SIMD
// --------------------------------------------------------

template <typename TScoreVec>
class Score<TScoreVec, ScoreSimdWrapper<Score<int32_t, PositionSpecificScoreSimd>>>
{
public:
    using TVecValue = typename Value<TScoreVec>::Type;
    using TBaseScoreSpec = typename Spec<Score<int32_t, PositionSpecificScoreSimd>>::Type;
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
    Score() = default;

    template <typename TScoreVal2, typename TScoreSpec2>
    explicit Score(Score<TScoreVal2, TScoreSpec2> const & pScore) :
        data_gap_extend(createVector<TScoreVec>(scoreGapExtend(pScore))),
        data_gap_open(createVector<TScoreVec>(scoreGapOpen(pScore))),
        _baseScore(pScore)
    {}
};

template <typename TValue, typename TVal1, typename TVal2>
TValue score(Score<TValue, ScoreSimdWrapper<Score<int32_t, PositionSpecificScoreSimd>>> const & sc,
             TVal1 const & val1,
             TVal2 const & val2)
{
    return score(sc._baseScore, val1, val2);
}

#endif

} // namespace seqan
