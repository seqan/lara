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

#include "data_types.hpp"
#include "io.hpp"
#include "lagrange.hpp"
#include "parameters.hpp"
#include "subgradient_solver.hpp"

int main (int argc, char const ** argv)
{
    // Parse arguments and options.
    lara::Parameters params(argc, argv);
    if (params.status != lara::Status::CONTINUE)
        return params.status == lara::Status::EXIT_OK ? 0 : 1;

    // Read input files and prepare structured sequences.
    lara::InputStorage store(params);
    size_t const problem_size = store.size() * (store.size() - 1) / 2;
    _VV(params, "Attempting to solve " << problem_size << " structural alignments.");
    _VV(params, store);
    lara::OutputTCoffeeLibrary tcLib(store);

    for (size_t idxA = 0ul; idxA < store.size() - 1ul; ++idxA)
    {
        for (size_t idxB = idxA + 1ul; idxB < store.size(); ++idxB)
        {
            _V(params, "SEQUENCE " << idxA << " WITH " << idxB);
            lara::Lagrange lagrange(store[idxA], store[idxB], params);
            lagrange.start();

            lara::SubgradientSolver solver(lagrange.getDimension(), params);
            lara::Status status = solver.solve(lagrange);
            if (status == lara::Status::EXIT_ERROR)
                return 1;

            tcLib.addAlignment(lagrange, idxA, idxB);
            if (problem_size == 1ul)
                lara::printAlignment(params.outFile, lagrange.getAlignment(), store[idxA].name, store[idxB].name);

        }
    }
    if (problem_size > 1ul)
        tcLib.print(params.outFile);
}
