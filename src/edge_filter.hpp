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

/*!\file edge_filter.hpp
 * \brief This file contains methods for filtering the alignment edges in the beginning of the LaRA algorithm.
 */

#include <algorithm>
#include <vector>

#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "data_types.hpp"
#include "score.hpp"

namespace lara
{

class PairwiseGotoh
{
private:
    size_t const lenA;
    size_t const lenB;

    std::vector<double> matrixM;
    std::vector<double> matrixH;
    std::vector<double> matrixV;

    //!\brief Shortcut for retrieving the specified matrix entry.
    double & get(std::vector<double> & matrix, size_t posA, size_t posB) const
    {
        return matrix[(lenB + 1ul) * posA + posB];
    }

public:
    explicit
    PairwiseGotoh(seqan::Rna5String const & seqA, seqan::Rna5String const & seqB, RnaScoreMatrix const & score)
        : lenA(seqan::length(seqA)), lenB(seqan::length(seqB))
    {
        matrixM.resize((lenA + 1) * (lenB + 1));
        matrixH.resize((lenA + 1) * (lenB + 1));
        matrixV.resize((lenA + 1) * (lenB + 1));
        double const go = score.data_gap_open;
        double const ge = score.data_gap_extend;

        // initialise DP matrix
        get(matrixM, 0, 0) = 0.0;
        get(matrixH, 0, 0) = negInfinity;
        get(matrixV, 0, 0) = negInfinity;


        for (size_t a = 0ul; a < lenA; ++a)
        {
            get(matrixM, a + 1, 0) = go + ge * a;
            get(matrixH, a + 1, 0) = negInfinity;
            get(matrixV, a + 1, 0) = go + ge * a;
        }

        for (size_t b = 0ul; b < lenB; ++b)
        {
            get(matrixM, 0, b + 1) = go + ge * b;
            get(matrixH, 0, b + 1) = go + ge * b;
            get(matrixV, 0, b + 1) = negInfinity;
        }

        // fill DP matrix
        for (size_t a = 0ul; a < lenA; ++a)
        {
            for (size_t b = 0ul; b < lenB; ++b)
            {
                get(matrixM, a + 1, b + 1) = std::max({get(matrixM, a, b),
                                                       get(matrixH, a, b),
                                                       get(matrixV, a, b)}) + seqan::score(score, seqA[a], seqB[b]);

                get(matrixH, a + 1, b + 1) = std::max({get(matrixM, a + 1, b) + go,
                                                       get(matrixH, a + 1, b) + ge,
                                                       get(matrixV, a + 1, b) + go});

                get(matrixV, a + 1, b + 1) = std::max({get(matrixM, a, b + 1) + go,
                                                       get(matrixH, a, b + 1) + go,
                                                       get(matrixV, a, b + 1) + ge});
            }
        }
    }

    double getPrefixScore(size_t posA, size_t posB)
    {
        assert(posA <= lenA && posB <= lenB);
        return std::max({get(matrixM, posA, posB), get(matrixH, posA, posB), get(matrixV, posA, posB)});
    }

    double getOptimalScore()
    {
        return getPrefixScore(lenA, lenB);
    }
};

void generateEdges(std::map<PosPair, size_t> & edges,  // OUT
                   seqan::Rna5String const & seqA,     // IN
                   seqan::Rna5String const & seqB,     // IN
                   RnaScoreMatrix const & scoreMatrix, // IN
                   double suboptimalDiff)              // IN
{
    seqan::ModifiedString<seqan::Rna5String const, seqan::ModReverse> reverseA(seqA);
    seqan::ModifiedString<seqan::Rna5String const, seqan::ModReverse> reverseB(seqB);
    PairwiseGotoh forward(seqA, seqB, scoreMatrix);
    PairwiseGotoh backward(reverseA, reverseB, scoreMatrix);
    assert(std::abs(forward.getOptimalScore() - backward.getOptimalScore()) < 1.e-3);
    double const threshold = forward.getOptimalScore() - suboptimalDiff;
    size_t const lenA = seqan::length(seqA);
    size_t const lenB = seqan::length(seqB);

    size_t edgeIdx = 0ul;
    for (size_t a = 0ul; a < lenA; ++a)
        for (size_t b = 0ul; b < lenB; ++b)
            if (forward.getPrefixScore(a, b)
                + seqan::score(scoreMatrix, seqA[a], seqB[b])
                + backward.getPrefixScore(lenA - a - 1, lenB - b - 1)
                >= threshold)
            {
                edges.emplace_hint(edges.end(), std::make_pair(a, b), edgeIdx++);
            }
}

} // namespace lara

