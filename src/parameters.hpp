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

/*!\file parameters.hpp
 * \brief This file contains all settings and parameters of LaRA.
 */

#include <string>
#include <thread>
#include <vector>

#include <seqan/arg_parse.h>
#include <seqan/score.h>

#include "data_types.hpp"
#include "rna_score_matrices.hpp"

namespace lara
{

template<typename TSequenceValue, typename TTag>
inline void setRnaScoreMatrix(seqan::Score<float, seqan::ScoreMatrix<TSequenceValue>> & matrix, TTag)
{
    float const * tab = RnaScoringMatrixData_<float, TSequenceValue, TTag>::getData();
    seqan::arrayCopy(tab, tab + RnaScoringMatrixData_<float, TSequenceValue, TTag>::TAB_SIZE, matrix.data_tab);
}

class Parameters
{
public:
    enum Status
    {
        EXIT_OK    = 0,
        EXIT_ERROR = 1,
        CONTINUE   = 2
    };
    Status status{};

    // GENERAL OPTIONS
    unsigned                 threads{};

    // INPUT OPTIONS
    std::string              inFile{};               // Name of input file
    std::string              refFile{};              // Name of input fileRef
    std::vector<std::string> dotplotFiles{};         // Names of dotplot files

    // OUTPUT OPTIONS
    std::string              outFile{};              // Name of output file (default: stdout)
    unsigned                 libraryScoreMin{};      // specify the minimum score for the T-Coffee library
    unsigned                 libraryScoreMax{};      // specify the maximum score for the T-Coffee library
    bool                     libraryScoreIsLinear{}; // whether T-Coffee scores are binary or linearly scaled

    // RUNTIME/QUALITY OPTIONS
    unsigned                 numIterations{};        // number of iterations
    unsigned                 maxNondecrIterations{}; // number of non-decreasing iterations
    float                    stepSizeFactor{};       // my, necessary for computing appropriate step sizes
    float                    epsilon{};              // max distance that means equality of upper and lower bound
    unsigned                 matching{};             // select matching algorithm
    float                    suboptimalDiff{};       // Gap open and extend costs for generating the alignment edges

    // SCORING OPTIONS
    float                    balance{};
    float                    sequenceScale{};        // scaling factor for the scores of the alignment edges
    unsigned                 structureScoring{};     // scoring mode for structures, either LOGARITHMIC or SCALE
    RnaScoreMatrix           rnaScore{};             // scoring matrix for scoring alignment edges (sequence score)

    // Constructor.
    Parameters(int argc, char const ** argv) noexcept
    {
        status = setParameters(argc, argv);
    }

private:
    inline Status setParameters(int argc, char const ** argv) noexcept
    {
        using namespace seqan;
        ArgumentParser parser;

        // Set up the parser.
        setAppName(parser, "lara");
        setShortDescription(parser, "Lagrangian Relaxed Alignment for RNA structures");
        setVersion(parser, "2.0.1");
        setDate(parser, "July 2019");
        addDescription(parser, "RNA structural alignment algorithm.");

        addUsageLine(parser, R"( -i \fIinFile\fP [\fIparameters\fP])");
        addUsageLine(parser, R"( -d \fIdpFile\fP -d \fIdpFile\fP [-d ...] [\fIparameters\fP])");

        addOption(parser, ArgParseOption("v", "verbose",
                                         "0: no additional outputs, 1: global statistics, "
                                         "2: extensive statistics for each batch of reads, 3: Debug output.",
                                         ArgParseArgument::INTEGER, "INT"));
        setMinValue(parser, "v", "0");
        setMaxValue(parser, "v", "3");
        setDefaultValue(parser, "v", "0");

        addOption(parser, ArgParseOption("j", "threads",
                                         "Use the number of specified threads. Decrease if memory problems occur. "
                                         "Value 0 tries to detect the maximum number.",
                                         ArgParseArgument::INTEGER, "INT"));
        setMinValue(parser, "j", "0");
        setDefaultValue(parser, "j", "1");

        // Input options
        addSection(parser, "Input Options");

        addOption(parser, ArgParseOption("i", "infile",
                                         "Path to the input file.",
                                         ArgParseArgument::INPUT_FILE, "IN"));

        addOption(parser, ArgParseOption("r", "reffile",
                                         "Path to the reference input file.",
                                         ArgParseArgument::INPUT_FILE, "IN"));

        addOption(parser, ArgParseOption("d", "dotplot",
                                         "Use dotplot files from RNAfold (*_dp.ps) as sequence and structure input.",
                                         ArgParseArgument::INPUT_FILE, "IN", true));

        // Output options
        addSection(parser, "Output Options");

        addOption(parser, ArgParseOption("w", "write",
                                         "Path to the output file. Default: stdout.",
                                         ArgParseArgument::OUTPUT_FILE, "OUT"));

        addOption(parser, ArgParseOption("l", "libscore",
                                         "The range of the scores for the T-Coffe library. "
                                         "Default: 500 1000 (binary).",
                                         ArgParseArgument::INTEGER, "MIN MAX",
                                         false, 2));

        // Runtime/Quality options
        addSection(parser, "Runtime/Quality Options");

        addOption(parser, ArgParseOption("n", "numiter",
                                         "The number of iterations.",
                                         ArgParseArgument::INTEGER, "INT"));
        setMinValue(parser, "n", "1");
        setDefaultValue(parser, "n", "500");

        addOption(parser, ArgParseOption("a", "maxnondecreasing",
                                         "The number of non-decreasing iterations.",
                                         ArgParseArgument::INTEGER, "INT"));
        setMinValue(parser, "a", "0");
        setDefaultValue(parser, "a", "50");

        addOption(parser, ArgParseOption("f", "factor",
                                         "The factor of how much to increase the step sizes.",
                                         ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(parser, "f", "1.0");

        addOption(parser, ArgParseOption("e", "epsilon",
                                         "Maximal distance that means equality of upper and lower bound.",
                                         ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(parser, "e", "0.01");

        addOption(parser, ArgParseOption("m", "matching",
                                         "Lookahead for greedy matching algorithm. Value 0 uses LEMON instead.",
                                         ArgParseArgument::INTEGER, "INT"));
        setMinValue(parser, "m", "0");
        setDefaultValue(parser, "m", "5");

        addOption(parser, ArgParseOption("u", "subopt",
                                         "Parameter for filtering alignment edges. Only those are created, whose prefix"
                                         " score + suffix score in the DP matrix is at most subopt below the "
                                         "optimal score.",
                                         ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(parser, "u", "40.0");


        // Scoring options
        addSection(parser, "Scoring Options");

        addOption(parser, ArgParseOption("b", "balance",
                                         "Factor of how much the sequence identity should influence the balance "
                                         "of sequence and structure score.",
                                         ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(parser, "b", "0.0");

        addOption(parser, ArgParseOption("c", "seqscale",
                                         "Scaling factor for the sequence scores (below 1 gives more impact for "
                                         "structure).",
                                         ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(parser, "c", "1.0");

        addOption(parser, ArgParseOption("p", "probscoremode",
                                         "The base pair probability scoring mode, either LOGARITHMIC (0), SCALE (1).",
                                         ArgParseArgument::INTEGER, "INT"));
        setMinValue(parser, "p", "0");
        setMaxValue(parser, "p", "1");
        setDefaultValue(parser, "p", "0");

        addOption(parser, ArgParseOption("x", "gapextend",
                                         "Gap extend costs for generating the alignment edges.",
                                         ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(parser, "x", "-2.0");

        addOption(parser, ArgParseOption("y", "gapopen",
                                         "Gap open costs for generating the alignment edges.",
                                         ArgParseArgument::DOUBLE, "FLOAT"));
        setDefaultValue(parser, "y", "-6.0");

        addOption(parser, ArgParseOption("s", "scorematrix",
                                         "The score matrix file for scoring alignment edges in the "
                                         "actual problem. Default: Ribosum65N.",
                                         ArgParseOption::STRING));

        // Parse.
        ArgumentParser::ParseResult parseResult = parse(parser, argc, argv);
        if (parseResult != ArgumentParser::ParseResult::PARSE_OK)
            return parseResult == ArgumentParser::ParseResult::PARSE_ERROR ? EXIT_ERROR : EXIT_OK;

        // Fill the parameters and check validity.

        // GENERAL OPTIONS
        getOptionValue(_VERBOSE_LEVEL, parser, "verbose");
        getOptionValue(threads, parser, "threads");
        if (threads == 0u)
        {
            unsigned nthreads = std::thread::hardware_concurrency();
            threads = nthreads != 0 ? nthreads : 1u;
        }

        // INPUT OPTIONS
        getOptionValue(inFile, parser, "infile");
        getOptionValue(refFile, parser, "reffile");
        dotplotFiles.resize(getOptionValueCount(parser, "dotplot"));
        for (size_t idx = 0ul; idx < dotplotFiles.size(); ++idx)
            getOptionValue(dotplotFiles[idx], parser, "dotplots", idx);

        if (inFile.empty() && dotplotFiles.empty())
        {
            printShortHelp(parser);
            return EXIT_ERROR;
        }

        // OUTPUT OPTIONS
        getOptionValue(outFile, parser, "write");
        if (isSet(parser, "write"))
        {
            _LOG(1, "The specified output file is " << outFile << std::endl);
        }

        getOptionValue(libraryScoreMin, parser, "libscore", 0);
        getOptionValue(libraryScoreMax, parser, "libscore", 1);
        libraryScoreIsLinear = isSet(parser, "libscore");

        // RUNTIME/QUALITY OPTIONS
        getOptionValue(numIterations, parser, "numiter");
        getOptionValue(maxNondecrIterations, parser, "maxnondecreasing");
        getOptionValue(stepSizeFactor, parser, "factor");
        getOptionValue(epsilon, parser, "epsilon");
        getOptionValue(matching, parser, "matching");
        getOptionValue(suboptimalDiff, parser, "subopt");

        // SCORING OPTIONS
        getOptionValue(balance, parser, "balance");
        getOptionValue(sequenceScale, parser, "seqscale");
        getOptionValue(structureScoring, parser, "probscoremode");
        getOptionValue(rnaScore.data_gap_open, parser, "gapopen");
        getOptionValue(rnaScore.data_gap_extend, parser, "gapextend");

        // set score matrix
        std::string scoreMatrixFile{};
        getOptionValue(scoreMatrixFile, parser, "scorematrix");
        if (scoreMatrixFile.empty())
        {
            _LOG(2, "Predefined Ribosum65 matrix will be used." << std::endl);
            setRnaScoreMatrix(rnaScore, Ribosum65N());
        }
        else if (loadScoreMatrix(rnaScore, scoreMatrixFile.c_str()))
        {
            _LOG(2, "Provided scoring matrix will be used: " << scoreMatrixFile << std::endl);
        }
        else
        {
            std::cerr << "Matrix file could not be opened: " << scoreMatrixFile
                      << "). Predefined Ribosum65 matrix will be used." << std::endl;
            setRnaScoreMatrix(rnaScore, Ribosum65N());
        }

        return Status::CONTINUE;
    }
};

} // namespace lara
