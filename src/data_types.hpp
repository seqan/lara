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

#include <iostream>
#include <limits>
#include <set>

#include <seqan/rna_io.h>
#include <seqan/score.h>

/*!
 * \brief Define verbosity levels.
 */
#define _V(_opt, _str) { if ((_opt).verbose > 0) std::cerr << _str << std::endl; }
#define _VV(_opt, _str) { if ((_opt).verbose > 1) std::cerr << _str << std::endl; }
#define _VVV(_opt, _str) { if ((_opt).verbose > 2) std::cerr << _str << std::endl; }

namespace lara
{

enum ScoringMode
{
    LOGARITHMIC,
    SCALE,
    ORIGINAL,
    RIBOSUM
};

enum TCoffeeMode
{
    PROPORTIONAL,
    SWITCH,
    ALLINTER,
    FIXEDINTER
};

enum Status
{
    EXIT_OK,
    EXIT_ERROR,
    CONTINUE
};

float const negInfinity = std::numeric_limits<float>::lowest();
float const posInfinity = std::numeric_limits<float>::max();

//! \brief Score Matrix type used in LaRA.
typedef seqan::Score<float, seqan::ScoreMatrix<seqan::Rna5>> RnaScoreMatrix;

//! \brief Pair of positions (usually in first and second sequence)
typedef std::pair<size_t, size_t>                           PosPair;
typedef std::pair<size_t, float>                            Contact;
typedef std::set<std::pair<float, size_t>>                  PriorityQueue;
typedef seqan::Align<seqan::Rna5String, seqan::ArrayGaps>   Alignment;
typedef seqan::Row<Alignment>::Type                         AlignmentRow;
typedef std::tuple<float, size_t, size_t>                   Interaction;    // probability, lineL, lineR
typedef std::set<Interaction>::iterator                     InteractionIterator;
typedef std::pair<InteractionIterator, InteractionIterator> EdgeConflict;
typedef std::set<InteractionIterator>                       InteractionSet;

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
