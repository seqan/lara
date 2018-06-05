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

/*!\file io.hpp
 * \brief This file contains LaRA's file handling, i.e. reading input files and writing output.
 */

#include <iostream>
#include <ostream>
#include <sstream>

#include <seqan/file.h>
#include <seqan/graph_types.h>
#include <seqan/rna_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "parameters.hpp"

extern "C" {
    #include <ViennaRNA/utils.h>
    #include <ViennaRNA/fold_vars.h>
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/part_func.h>
}

namespace lara
{

class InputStorage : public std::vector<seqan::RnaRecord>
{
public:
    explicit InputStorage(Parameters const & params)
    {
        readRnaFile(params.inFile);
        readRnaFile(params.inFileRef);
        _VV(params, "Successfully read " << size() << " records.");

        if (params.dotplotFile.size() == size())
        {
            // Load base pair probabilities from dot plot file.
            for (size_t fileIdx = 0u; fileIdx < size(); ++fileIdx)
                extractBppFromDotplot(at(fileIdx), params.dotplotFile[fileIdx]);
            _VV(params, "Successfully extracted base pair probabilities from given dotplot files.");
        }
        else
        {
            // If not present, compute the weighted interaction edges using ViennaRNA functions.
            bool const logScoring = params.structureScoring == ScoringMode::LOGARITHMIC;
            bool usedVienna = false;
            for (seqan::RnaRecord & record : *this)
                computeStructure(record, usedVienna, logScoring);
            if (usedVienna)
                _VV(params, "Computed missing base pair probabilities with ViennaRNA library.");
        }
    }

private:
    void readRnaFile(seqan::CharString filename)
    {
        if (seqan::empty(filename))
            return;

        seqan::RnaStructFileIn rnaStructFile;
        if (seqan::open(rnaStructFile, toCString(filename), seqan::OPEN_RDONLY))
        {
            seqan::RnaHeader header;                   // dummy, just for Ebpseq files
            seqan::readHeader(header, rnaStructFile);
            seqan::RnaRecord rec;
            while (!seqan::atEnd(rnaStructFile))
            {
                seqan::readRecord(rec, rnaStructFile);
                push_back(rec);
            }
            seqan::close(rnaStructFile);
        }
        else
        {
            // Read the file.
            seqan::SeqFileIn seqFileIn(toCString(filename));
            seqan::StringSet<seqan::CharString>  ids;
            seqan::StringSet<seqan::IupacString> seqs;
            seqan::StringSet<seqan::CharString>  quals;
            seqan::readRecords(ids, seqs, quals, seqFileIn);
            seqan::close(seqFileIn);

            // Fill the data structures: identifier and sequence.
            reserve(size() + seqan::length(ids));
            SEQAN_ASSERT_EQ(seqan::length(ids), seqan::length(seqs));
            for (size_t idx = 0ul; idx < seqan::length(ids); ++idx)
            {
                seqan::RnaRecord rec{};
                rec.name     = ids[idx];
                rec.sequence = seqan::convert<seqan::Rna5String>(seqs[idx]);

                // For FastQ files: add quality annotation.
                if (seqan::length(quals) == seqan::length(ids))
                    rec.quality = quals[idx];
                push_back(rec);
            }
        }
    }

    void extractBppFromDotplot(seqan::RnaRecord & rnaRecord, std::string const & dotplotFile)
    {
        using namespace seqan;

        double const minProb = 0.003; // taken from LISA > Lara

        // add vertices to graph
        RnaStructureGraph bppMatrGraph;
        for (size_t idx = 0u; idx < length(rnaRecord.sequence); ++idx)
            addVertex(bppMatrGraph.inter);

        // open dotplot file and read lines
        std::ifstream file(dotplotFile);
        std::string   line;
        while (std::getline(file, line))
        {
            if (line.find("ubox") == std::string::npos)
                continue;

            std::istringstream iss(line);
            unsigned           iPos, jPos;
            double             prob;
            if (iss >> iPos >> jPos >> prob) // read values from line
            {   // create edges for graph
                SEQAN_ASSERT(iPos > 0 && iPos <= length(rnaRecord.sequence));
                SEQAN_ASSERT(jPos > 0 && jPos <= length(rnaRecord.sequence));
                // convert indices from range 1..length to 0..length-1
                if (prob * prob > minProb) // dot plot contains sqrt(prob)
                    addEdge(bppMatrGraph.inter, iPos - 1, jPos - 1, log(prob * prob / minProb));
            }
        }
        bppMatrGraph.specs = CharString("ViennaRNA dot plot from file " + std::string(dotplotFile));
        append(rnaRecord.bppMatrGraphs, bppMatrGraph);
    }

    void computeStructure(seqan::RnaRecord & rnaRecord, bool & usedVienna, bool logStructureScoring)
    {
        if (!seqan::empty(rnaRecord.bppMatrGraphs))
            return;

        usedVienna = true;
        size_t const length = seqan::length(rnaRecord.sequence);
        seqan::String<char, seqan::CStyle> sequence{rnaRecord.sequence};

        // Compute the partition function and base pair probabilities with ViennaRNA.
        seqan::RnaStructureGraph bppMatrGraph;
        init_pf_fold(static_cast<int>(length));
        bppMatrGraph.energy = pf_fold(seqan::toCString(sequence), nullptr);
        bppMatrGraph.specs = seqan::CharString{"ViennaRNA pf_fold"};

        for (size_t idx = 0u; idx < length; ++idx)
            seqan::addVertex(bppMatrGraph.inter);

        double const  minProb = 0.003; // taken from LISA > Lara
        for (size_t i = 0u; i < length; ++i)
        {
            for (size_t j = i + 1u; j < length; ++j)
            {
                if (logStructureScoring)
                {
                    if (pr[iindx[i + 1] - (j + 1)] > minProb)
                        seqan::addEdge(bppMatrGraph.inter, i, j, log(pr[iindx[i + 1] - (j + 1)] / minProb));
                }
                else
                {
                    if (pr[iindx[i + 1] - (j + 1)] > 0.0)
                        seqan::addEdge(bppMatrGraph.inter, i, j, pr[iindx[i + 1] - (j + 1)]);
                }
            }
        }
        seqan::append(rnaRecord.bppMatrGraphs, bppMatrGraph);

        // Compute the fixed structure with ViennaRNA.
        auto * structure = new char[length + 1];
        initialize_fold(static_cast<int>(length));
        float energy = fold(seqan::toCString(sequence), structure);
        seqan::bracket2graph(rnaRecord.fixedGraphs, seqan::CharString{structure}); // appends the graph
        seqan::back(rnaRecord.fixedGraphs).energy = energy;
        seqan::back(rnaRecord.fixedGraphs).specs = seqan::CharString{"ViennaRNA fold"};
        delete[] structure;
    }
};

std::ostream & operator<<(std::ostream & stream, InputStorage & store)
{
    using namespace seqan;
    for (seqan::RnaRecord const & rec : store)
    {
        writeRecord(stream, rec, DotBracket());
    }
    return stream;
}

} // namespace lara

