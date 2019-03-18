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
#include <ostream>
#include <sstream>

#include <seqan/consensus.h>
#include <seqan/score.h>

namespace seqan
{

struct RnaStructureScore;

template <>
class Score<int32_t, RnaStructureScore>
{
public:
    std::vector<std::vector<int32_t>> matrix;
    int32_t data_gap_open;
    int32_t data_gap_extend;
};

template <typename TSequence>
struct SequenceEntryForScore<Score<int32_t, RnaStructureScore>, TSequence>
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template <typename TSeq1, typename TSeq2>
inline
int32_t score(Score<int32_t, RnaStructureScore> const & sc,
              ConsensusScoreSequenceEntry<TSeq1> const & entryH,
              ConsensusScoreSequenceEntry<TSeq2> const & entryV)
{
    return sc.matrix[position(entryH)][position(entryV)];
}

template <typename TSequence, typename TPosition>
inline
ConsensusScoreSequenceEntry<TSequence> sequenceEntryForScore(Score<int32_t, RnaStructureScore> const &,
                                                             TSequence const & seq,
                                                             TPosition pos)
{
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

} // namespace seqan

