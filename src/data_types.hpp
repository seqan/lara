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

/*!\file data_types.hpp
 * \brief This file contains data structures and macros for LaRA.
 */

#include <chrono>
#include <iostream>
#include <functional>
#include <limits>
#include <set>

#include <seqan/rna_io.h>
#include <seqan/score.h>

#define _LOG(_level, _str)   { if (lara::_VERBOSE_LEVEL >= (_level)) std::cerr << _str; }

namespace lara
{

/*!
 * \brief Define verbosity levels.
 */
int _VERBOSE_LEVEL = 0;

enum ScoringMode
{
    LOGARITHMIC,
    SCALE
};

//enum TCoffeeMode
//{
//    PROPORTIONAL,
//    SWITCH,
//    ALLINTER,
//    FIXEDINTER
//};

//! \brief Pair of positions (usually in first and second sequence)
typedef int32_t                                                  ScoreType;
typedef uint32_t                                                 UnsignedType;
typedef seqan::Score<ScoreType, seqan::ScoreMatrix<seqan::Rna5>> SeqScoreMatrix;
typedef std::pair<size_t, size_t>                                PosPair;
typedef std::pair<ScoreType, size_t>                             Contact;
typedef std::set<Contact>                                        PriorityQueue;
typedef seqan::Gaps<seqan::String<unsigned>, seqan::ArrayGaps>   GappedSeq;
typedef std::pair<GappedSeq, GappedSeq>                          Alignment;
typedef std::pair<PosPair, std::vector<std::tuple<size_t, size_t, unsigned>>> WeightedAlignedColumns;
typedef std::chrono::steady_clock                                Clock;

ScoreType const infinity = std::numeric_limits<ScoreType>::max() / 3 * 2;
float const factor2int = 8192.f;

// Time helper functions.
template <typename duration_unit = std::chrono::milliseconds>
inline typename duration_unit::rep timeDiff(Clock::time_point start)
{
    return std::chrono::duration_cast<duration_unit>(Clock::now() - start).count();
}

} // namespace lara
