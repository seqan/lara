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

#include <algorithm>
#include <cctype>
#include <iostream>
#include <ostream>
#include <sstream>

#include <seqan/file.h>
#include <seqan/graph_types.h>
#include <seqan/rna_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "lagrange.hpp"
#include "parameters.hpp"

#ifdef VIENNA_RNA_FOUND
extern "C" {
    #include <ViennaRNA/utils.h>
    #include <ViennaRNA/fold_vars.h>
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/part_func.h>
}
#endif

namespace lara
{

class InputStorage : public std::vector<seqan::RnaRecord>
{
public:
    explicit InputStorage(Parameters const & params)
    {
        readRnaFile(params.inFile);
        readRnaFile(params.inFileRef);

        // If not present, compute the weighted interaction edges using ViennaRNA functions.
        bool const logScoring = params.structureScoring == ScoringMode::LOGARITHMIC;
        bool usedVienna = false;
        for (seqan::RnaRecord & record : *this)
            computeStructure(record, usedVienna, logScoring);
        if (usedVienna)
            _VV(params, "Computed missing base pair probabilities with ViennaRNA library.");

        if (!params.dotplotFile.empty())
        {
            // Load base pair probabilities from dot plot file.
            for (std::string const & filename : params.dotplotFile)
            {
                seqan::RnaRecord rec;
                extractBppFromDotplot(rec, filename);
                push_back(rec);
            }
            _VV(params, "Successfully extracted base pair probabilities from given dotplot files.");
        }

        if (size() <= 1)
            throw std::runtime_error("ERROR: The given file(s) must contain at least two sequences.");

        for (seqan::RnaRecord & record : *this)
            _VVV(params, record.bppMatrGraphs[0].inter);
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

    void extractBppFromDotplot(seqan::RnaRecord & rnaRecord, std::string const & filename)
    {
        using namespace seqan;

        double const minProb = 0.003; // taken from LISA > Lara

        // open dotplot file and read lines
        std::ifstream file(filename);
        std::string   line;
        unsigned      iPos, jPos;
        double        prob;

        if (!file.is_open())
            throw std::runtime_error("ERROR: Cannot open file " + filename);

        while (std::getline(file, line))
        {
            if (line.find("/sequence") != std::string::npos)
            {
                while (std::getline(file, line) && line.find(')') == std::string::npos)
                {
                    line.erase(std::remove_if(line.begin(), line.end(),
                                              [] (unsigned char x) { return std::isalpha(x) == 0; }),
                               line.end());
                    seqan::append(rnaRecord.sequence, line);
                    std::cerr << line << std::endl;
                }
                break;
            }
        }
        SEQAN_ASSERT(!seqan::empty(rnaRecord.sequence));

        // add vertices to graph
        RnaStructureGraph bppMatrGraph;
        RnaStructureGraph fixedGraph;
        for (size_t idx = 0u; idx < length(rnaRecord.sequence); ++idx)
        {
            addVertex(bppMatrGraph.inter);
            addVertex(fixedGraph.inter);
        }

        while (std::getline(file, line))
        {
            if (line.find("ubox") != std::string::npos)
            {
                std::istringstream iss(line);
                if (iss >> iPos >> jPos >> prob) // read values from line
                {   // create edges for graph
                    SEQAN_ASSERT(iPos > 0 && iPos <= length(rnaRecord.sequence));
                    SEQAN_ASSERT(jPos > 0 && jPos <= length(rnaRecord.sequence));
                    // convert indices from range 1..length to 0..length-1
                    if (prob * prob > minProb) // dot plot contains sqrt(prob)
                        addEdge(bppMatrGraph.inter, iPos - 1, jPos - 1, log(prob * prob / minProb));
                }
            }
            else if (line.find("lbox") != std::string::npos)
            {
                std::istringstream iss(line);
                if (iss >> iPos >> jPos >> prob) // read values from line
                {   // create edges for graph
                    SEQAN_ASSERT(iPos > 0 && iPos <= length(rnaRecord.sequence));
                    SEQAN_ASSERT(jPos > 0 && jPos <= length(rnaRecord.sequence));
                    // convert indices from range 1..length to 0..length-1
                    addEdge(fixedGraph.inter, iPos - 1, jPos - 1, 1.0);
                }
            }
        }
        if (seqan::numEdges(bppMatrGraph.inter) == 0)
            throw std::runtime_error("WARNING: No structure information found in file " + filename);

        bppMatrGraph.specs = CharString("ViennaRNA dot plot from file " + std::string(filename));
        fixedGraph.specs   = CharString("ViennaRNA dot plot from file " + std::string(filename));

        std::string name = filename.substr(filename.find_last_of("/\\") + 1);
        rnaRecord.name = name.substr(0, name.rfind(".ps")).substr(0, name.rfind("_dp"));
        append(rnaRecord.bppMatrGraphs, bppMatrGraph);
        append(rnaRecord.fixedGraphs, fixedGraph);
    }

    void computeStructure(seqan::RnaRecord & rnaRecord, bool & usedVienna, bool logStructureScoring)
    {
        if (!seqan::empty(rnaRecord.bppMatrGraphs))
            return;

#ifndef VIENNA_RNA_FOUND
        std::cerr << "Cannot compute a structure without the ViennaRNA library. Please install ViennaRNA and try again."
                  << std::endl;
        exit(1);
#endif

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

// FASTA output

inline
void printAlignment(std::ostream & stream, Alignment const & alignment, seqan::CharString const & nameA, seqan::CharString const & nameB)
{
    stream << ">" << nameA << std::endl << seqan::row(alignment, 0) << std::endl;
    stream << ">" << nameB << std::endl << seqan::row(alignment, 1) << std::endl;
}

inline
void printAlignment(seqan::CharString & filename, Alignment const & alignment, seqan::CharString const & nameA, seqan::CharString const & nameB)
{
    if (seqan::empty(filename))
    {
        printAlignment(std::cout, alignment, nameA, nameB);
    }
    else
    {
        std::ofstream fastaFile;
        fastaFile.open(seqan::toCString(filename), std::ios::out);
        if (fastaFile.is_open())
        {
            printAlignment(fastaFile, alignment, nameA, nameB);
            fastaFile.close();
        }
        else
        {
            std::cerr << "Unable to open the specified output file for writing: " << filename << std::endl;
        }
    }
}

// T-COFFEE output

class OutputTCoffeeLibrary
{
private:
    std::ostringstream sstream;
    size_t const numSequences;

public:
    explicit OutputTCoffeeLibrary(InputStorage const & data) : numSequences(data.size())
    {
        if (numSequences > 2)
        {
            sstream << "! T-COFFEE_LIB_FORMAT_01" << std::endl;
            sstream << numSequences << std::endl;
            for (seqan::RnaRecord const & rec : data)
            {
                sstream << rec.name << " " << seqan::length(rec.sequence) << " " << rec.sequence << std::endl;
            }
        }
    }

    void addAlignment(Lagrange const & lagrange, size_t seqIndexA, size_t seqIndexB)
    {
        if (numSequences > 2)
        {
            sstream << "# " << (seqIndexA + 1) << " " << (seqIndexB + 1) << std::endl;
            for (auto const & elem : lagrange.getStructureLines())
            {
                sstream << std::get<0>(elem) << " " << std::get<1>(elem) << " " << (std::get<2>(elem) ? 1000 : 500)
                        << std::endl;
            }
        }
    }

    friend std::ostream & operator<<(std::ostream & stream, OutputTCoffeeLibrary & library);

    void print(seqan::CharString & filename)
    {
        if (seqan::empty(filename))
        {
            std::cout << *this;
        }
        else
        {
            std::ofstream tcLibFile;
            tcLibFile.open(seqan::toCString(filename), std::ios::out);
            if (tcLibFile.is_open())
            {
                tcLibFile << *this;
                tcLibFile.close();
            }
            else
            {
                std::cerr << "Unable to open the specified output file for writing: " << filename << std::endl;
            }
        }
    }

};

std::ostream & operator<<(std::ostream & stream, OutputTCoffeeLibrary & library)
{
    return stream << library.sstream.str() << "! SEQ_1_TO_N" << std::endl;
}

} // namespace lara

