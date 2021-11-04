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

/*!\file io.hpp
 * \brief This file contains LaRA's file handling, i.e. reading input files and writing output.
 */

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <vector>

#include <seqan/file.h>
#include <seqan/graph_types.h>
#include <seqan/rna_io.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "data_types.hpp"
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
    explicit InputStorage(Parameters const & params) : err(false)
    {
        _LOG(1, "2) Read input files...\n");
        Clock::time_point timeRead = Clock::now();
        readRnaFile(params.inFile);
        readRnaFile(params.refFile);
        if (err)
            return;
        _LOG(1, "   * sequence/structure files -> " << timeDiff(timeRead) << "ms\n");

        // If not present, compute the weighted interaction edges using ViennaRNA functions.
        Clock::time_point timeRnaFold = Clock::now();
        bool const logScoring = params.structureScoring == ScoringMode::LOGARITHMIC;
        bool usedVienna = false;
        for (seqan::RnaRecord & record : *this)
            computeStructure(record, usedVienna, logScoring);
        if (usedVienna)
            _LOG(1, "   * compute missing base pair probabilities with RNAfold -> " << timeDiff(timeRnaFold) << "ms\n");

        if (!params.dotplotFiles.empty())
        {
            // Load base pair probabilities from dot plot file.
            Clock::time_point timeDotplot = Clock::now();
            for (std::string const & filename : params.dotplotFiles)
            {
                seqan::RnaRecord rec;
                err = extractBppFromDotplot(rec, filename);
                if (!err && seqan::empty(rec.bppMatrGraphs))
                {
                    std::cerr << "ERROR: The dotplot file " << filename << " does not contain any base pair "
                                 "probabilities. Please make sure that you execute RNAfold with -p option and specify "
                                 "the resulting _dp.ps file!\n";
                    err = true;
                }
                if (err)
                    return;
                push_back(rec);
            }
            _LOG(1, "   * dotplot files -> " << timeDiff(timeDotplot) << "ms\n");
        }

        if (size() <= 1)
        {
            std::cerr << "ERROR: The given file(s) must contain at least two sequences.\n";
            err = true;
        }
        else if (size() > 2 && params.outFormat == "fasta")
        {
            std::cerr << "WARNING: We are computing more than one pairwise alignment "
                         "and you have selected fasta output.\n";
        }
    }

    bool had_err() const
    {
        return err;
    }

private:
    bool err;

    void readRnaFile(std::string const & filename)
    {
        if (filename.empty())
            return;

        seqan::RnaStructFileIn rnaStructFile;
        seqan::SeqFileIn seqFileIn;
        if (seqan::open(rnaStructFile, filename.c_str(), seqan::OPEN_RDONLY))
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
        else if (seqan::open(seqFileIn, filename.c_str(), seqan::OPEN_RDONLY))
        {
            // Read the file.
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
        else // try to read Fasta
        {
            uint32_t recCount = 0;
            std::ifstream input(filename);
            if (!input.is_open())
            {
                std::cerr << "ERROR: The file " + filename + " cannot be opened.\n";
                err = true;
                return;
            }
            std::istreambuf_iterator<char> iter(input);
            std::istreambuf_iterator<char> end_iter;
            while (!iter.equal(end_iter))
            {
                seqan::RnaRecord rec{};
                seqan::IupacString sequence;
                seqan::readRecord(rec.name, sequence, iter, seqan::Fasta());
                rec.recordID = recCount++;
                rec.sequence = seqan::convert<seqan::Rna5String>(sequence);
                push_back(rec);
            }
        }
    }

    static bool extractBppFromDotplot(seqan::RnaRecord & rnaRecord, std::string const & filename)
    {
        using namespace seqan;

        float const minProb = 0.003f; // taken from LISA > Lara

        // open dotplot file and read lines
        std::ifstream file(filename);
        std::string line;
        unsigned iPos{};
        unsigned jPos{};
        float prob{};

        if (!file.is_open())
        {
            std::cerr << "ERROR: Cannot open dotplot file " << filename << std::endl;
            return true;
        }

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

        bppMatrGraph.specs = CharString{"ViennaRNA dot plot from file " + filename};
        fixedGraph.specs   = CharString{"ViennaRNA structure from file " + filename};

        std::string name = filename.substr(filename.find_last_of("/\\") + 1);
        rnaRecord.name = name.substr(0, name.rfind(".ps")).substr(0, name.rfind("_dp"));
        if (seqan::numEdges(bppMatrGraph.inter) > 0)
            append(rnaRecord.bppMatrGraphs, bppMatrGraph);
        else
            append(rnaRecord.fixedGraphs, fixedGraph);
        return false;
    }

    static void computeStructure(seqan::RnaRecord & rnaRecord, bool & usedVienna, bool logStructureScoring)
    {
        if (!seqan::empty(rnaRecord.bppMatrGraphs))
            return;
        if (!seqan::empty(rnaRecord.fixedGraphs))
        {
            // we only have a fixed graph: increase the importance of the edges
            typedef typename seqan::Iterator<seqan::Graph<seqan::Undirected<double>>, seqan::EdgeIterator>::Type EdgeIt;
            for (EdgeIt edgeIt(seqan::front(rnaRecord.fixedGraphs).inter); !seqan::atEnd(edgeIt); seqan::goNext(edgeIt))
                seqan::cargo(*edgeIt) *= 10;
            return;
        }

#ifdef VIENNA_RNA_FOUND
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

        float const  minProb = 0.003f; // taken from LISA > Lara
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
                    if (pr[iindx[i + 1] - (j + 1)] > minProb)
                        seqan::addEdge(bppMatrGraph.inter, i, j, pr[iindx[i + 1] - (j + 1)]);
                }
            }
        }
        seqan::append(rnaRecord.bppMatrGraphs, bppMatrGraph);
#else
        std::cerr << "Cannot compute a structure without the ViennaRNA library. Please install ViennaRNA and try again."
                  << std::endl;
        (void) usedVienna;
        (void) logStructureScoring;
        exit(1);
#endif
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

// sort alignment order for the first sequence length
struct CompareSeqLength
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

// T-COFFEE or pairwise gapped sequence output
class OutputLibrary
{
private:
    InputStorage const & data;
    std::set<std::pair<WeightedAlignedColumns, ScoreType>> alignments{};
    std::string const format;

public:
    explicit OutputLibrary(InputStorage const & input, std::string outputFormat) :
        data(input), format(std::move(outputFormat))
    {}

    void addAlignment(WeightedAlignedColumns const & structureLines, ScoreType const alignmentScore)
    {
        alignments.emplace(structureLines, alignmentScore);
    }

    friend std::ostream & operator<<(std::ostream & stream, OutputLibrary const & library);

    void printLib(std::ostream & stream) const
    {
        stream << "! T-COFFEE_LIB_FORMAT_01\n" << data.size() << '\n';
        for (seqan::RnaRecord const & rec : data)
        {
            stream << rec.name << " " << seqan::length(rec.sequence) << " " << rec.sequence << '\n';
        }

        for (auto const & ali : alignments)
        {
            WeightedAlignedColumns const & structureLines = ali.first;
            stream << "# " << (structureLines.first.first + 1) << " " << (structureLines.first.second + 1) << '\n';
            for (std::tuple<size_t, size_t, unsigned> const & elem : structureLines.second)
                stream << std::get<0>(elem) + 1 << " " << std::get<1>(elem) + 1 << " " << std::get<2>(elem) << '\n';
        }

        stream << "! SEQ_1_TO_N\n";
    }

    void printAlignments(std::ostream & stream) const
    {
        for (auto const & ali : alignments)
        {
            WeightedAlignedColumns const & structureLines = ali.first;
            seqan::RnaRecord const & rec1 = data[structureLines.first.first];
            seqan::RnaRecord const & rec2 = data[structureLines.first.second];
            std::pair<std::ostringstream, std::ostringstream> gapped{};
            PosPair curr{0, 0};

            for (auto const & column : structureLines.second)
            {
                while (curr.first < std::get<0>(column))
                {
                    gapped.first << rec1.sequence[curr.first++];
                    gapped.second << '-';
                }

                while (curr.second < std::get<1>(column))
                {
                    gapped.first << '-';
                    gapped.second << rec2.sequence[curr.second++];
                }

                gapped.first << rec1.sequence[curr.first++];
                gapped.second << rec2.sequence[curr.second++];
            }
            while (curr.first < seqan::length(rec1.sequence))
            {
                gapped.first << rec1.sequence[curr.first++];
                gapped.second << '-';
            }
            while (curr.second < seqan::length(rec2.sequence))
            {
                gapped.first << '-';
                gapped.second << rec2.sequence[curr.second++];
            }

            if (format == "pairs")
            {
                stream << '>' << rec1.name << " && " << rec2.name << " (score " << ali.second / factor2int << ")\n";
                stream << gapped.first.str() << '\n' << gapped.second.str() << '\n';
            }
            else // format == "fasta"
            {
                stream << '>' << rec1.name << '\n' << gapped.first.str() << '\n';
                stream << '>' << rec2.name << '\n' << gapped.second.str() << '\n';
            }
        }
    }

    void print(std::ostream & stream) const
    {
        if (format == "lib")
            printLib(stream);
        else
            printAlignments(stream);
    }

    void print(std::string const & filename) const
    {
        _LOG(1, "4) Write results...\n");
        Clock::time_point timePrint = Clock::now();
        if (filename.empty())
        {
            print(std::cout);
            _LOG(1, "   * to stdout -> " << timeDiff(timePrint) << "ms\n");
        }
        else
        {
            std::ofstream file;
            file.open(filename.c_str(), std::ios::out);
            if (file.is_open())
            {
                print(file);
                file.close();
                _LOG(1, "   * to file " << filename << " -> " << timeDiff(timePrint) << "ms\n");
            }
            else
            {
                std::cerr << "Error: Unable to open the file for writing: " << filename << '\n';
            }
        }
    }
};

std::ostream & operator<<(std::ostream & stream, OutputLibrary const & library)
{
    library.print(stream);
    return stream;
}

} // namespace lara
