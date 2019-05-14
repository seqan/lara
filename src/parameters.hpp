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

#include <thread>

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
    seqan::arrayCopy(tab, tab + RnaScoringMatrixData_<float, TSequenceValue, TTag>::TAB_SIZE,
                     matrix.data_tab);
}

class Parameters
{
public:
    // Name of input file
    seqan::CharString        inFile{};
    // Name of input fileRef
    seqan::CharString        inFileRef{};
    // Name of dotplot file
    std::vector<std::string> dotplotFile{};
    // Name of output file (default: stdout)
    seqan::CharString        outFile{};
    // temporary directory where to save intermediate files. Default: use the input file directory.
    seqan::CharString        tmpDir{};
    // Use the gap scheme to be used in the N-W alignment (default: affine(0))
    unsigned                 affineLinearDgs{0u};
    // Use the global local or global-Unconstrained algorithm (default: global(0) - local(1) )
    bool                     alignLocally{false};
    // Parameter used during the RNAfold execution to select the minimum energy to be considered
    float                    thrBppm{1e-15f}; // old Lara: 0.1
    // number of iterations
    unsigned                 numIterations{500u};
    // number of non-decreasing iterations
    unsigned                 maxNondecrIterations{50u};
    // value to be considered for the equality of upper and lower bounds difference
    float                    epsilon{0.01f};
    // my, necessary for computing appropriate step sizes
    float                    stepSizeFactor{1.0f};
    // scoring matrix name that should be used for scoring alignment edges in the actual problem
    seqan::CharString        laraScoreMatrixName{};
    //    Score<float, ScoreMatrix<Rna5, Default> > laraScoreMatrix;
    RnaScoreMatrix           laraScoreMatrix;
    //    TScoringSchemeRib laraScoreMatrixRib;
    // Gap open and extend costs for generating the alignment edges
    float                    suboptimalDiff{40.0f};
    // Gap open and extend costs for generating the alignment edges
    float                    laraGapOpen{-6.0f};
    float                    laraGapExtend{-2.0f};
    // scaling factor for the scores of the alignment edges
    // Specifies the contribution of the sequence scores (specified by the larascore matrix) to the overall structural
    // alignment.
    float                    sequenceScale{1.0f};
    float                    balance{0.0f};
    // gap penalty for RSA
    //    float rsaGapPenalty{3.0};
    // scoring mode, either LOGARITHMIC, SCALE, ORIGINAL, RIBOSUM
    unsigned                 structureScoring{ScoringMode::LOGARITHMIC};
    // define the weight of _half_ an interaction match for fixed structures
    float                    fixedStructWeight{8.0f};
    // if structureScoring=SCALING then we have to give a scaling factor
    float                    scalingFactor{1.0f};
    // specify the location of T-COFFEE
    seqan::CharString        tcoffeeLocation{"t_coffee/t_coffee_5.05"};
    // specify the method to be used to create the T-Coffe library
    unsigned                 tcoffeeLibMode{0};
    Status                   status;
    unsigned                 matching{5u};
    unsigned                 num_threads{1u};

    Parameters(int argc, char const ** argv)
    {
        status = setParameters(argc, argv);
    }

private:
    inline Status setParameters(int argc, char const ** argv)
    {
        using namespace seqan;
        ArgumentParser parser;

        setAppName(parser, "lara");
        setShortDescription(parser, "Lagrangian Relaxed Alignment for RNA structures");
        setVersion(parser, "2.0.0");
        setDate(parser, "January 2019");
        addDescription(parser, "RNA structural alignment algorithm.");

        addUsageLine(parser, " -i \\fIinFile\\fP [\\fIparameters\\fP]");
        addUsageLine(parser, " -d \\fIdpFile\\fP -d \\fIdpFile\\fP [-d ...] [\\fIparameters\\fP]");

        addOption(parser, ArgParseOption("v", "verbose",
                                         "0: no additional outputs, 1: global statistics, "
                                         "2: extensive statistics for each batch of reads, 3: Debug output. (0)",
                                         ArgParseArgument::INTEGER, "INT"));

        // Input options
        addSection(parser, "Input Options");

        addOption(parser, ArgParseOption("i", "inFile",
                                         "Path to the input file",
                                         ArgParseArgument::INPUT_FILE, "IN"));

        addOption(parser, ArgParseOption("ir", "inFileRef",
                                         "Path to the reference input file",
                                         ArgParseArgument::INPUT_FILE, "IN"));

        addOption(parser, ArgParseOption("d", "dotplots",
                                         "Use dotplot files from RNAfold (*_dp.ps) as sequence and structure input.",
                                         ArgParseArgument::INPUT_FILE, "IN", true));

        // Output options
        addSection(parser, "Output Options");

        addOption(parser, ArgParseOption("w", "outFile",
                                         "Path to the output file (default: stdout)",
                                         ArgParseArgument::OUTPUT_FILE, "OUT"));

        addOption(parser, ArgParseOption("td", "tmpDir",
                                         "A temporary directory where to save intermediate files. (input file directory)",
                                         ArgParseOption::STRING));

        addOption(parser, ArgParseOption("tcl", "tcoffeeLocation",
                                         "location of T-COFFEE.",
                                         ArgParseOption::STRING));

        addOption(parser, ArgParseOption("tcm", "tcoffeeLibMode",
                                         "Method used to score the T-Coffe library. Either 0: switch(500/1000) or "
                                         "NUM: proportional in [NUM-250..NUM+250]. (0)",
                                         ArgParseArgument::INTEGER, "INT"));

        // Alignment options
        addSection(parser, "LaRA Alignment Options");

        addOption(parser, ArgParseOption("g", "affineLinearDgs",
                                         "Chose the gap scheme affine(0) linear(1) or dynamic(2) "
                                         "to be used in the alignment. (affine(0)). ",
                                         ArgParseArgument::INTEGER, "INT"));

        addOption(parser, ArgParseOption("a", "local",
                                         "Perform local alignment. (False)"));

        addOption(parser, ArgParseOption("tb", "thrBppm",
                                         "(Parameter used during the RNAfold execution to select the minimum energy to be "
                                         "considered (1e-15)",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("iter", "iterations",
                                         "number of iterations. ",
                                         ArgParseArgument::INTEGER, "INT"));

        addOption(parser, ArgParseOption("nditer", "maxNondecreasingIterations",
                                         "number of non-decreasing iterations. (50)",
                                         ArgParseArgument::INTEGER, "INT"));

        addOption(parser, ArgParseOption("ep", "epsilon",
                                         "value to be considered for the equality of upper and lower bounds difference",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("my", "stepSizeFactor",
                                         "necessary for computing appropriate step sizes.",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("lsm", "laraScoreMatrixName",
                                         "scoring matrix name that should be used for scoring alignment edges in the "
                                         "actual problem",
                                         ArgParseOption::STRING));

        addOption(parser, ArgParseOption("s", "suboptimalDiff",
                                         "Parameter for filtering alignment edges. Only those are created, whose prefix"
                                         " score + suffix score in the DP matrix is at most suboptimalDiff below the "
                                         "optimal score.",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("lgo", "laraGapOpen",
                                         "Gap open costs for generating the alignment edges",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("lge", "laraGapExtend",
                                         "Gap extend costs for generating the alignment edges",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("ssc", "sequenceScale",
                                         "Scaling factor for the sequence scores (below 1 gives more impact for "
                                         "structure).",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("b", "balance",
                                         "Factor of how much the sequence identity should influence the balance "
                                         "of sequence and structure score. (0)",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("stsc", "structureScoring",
                                         "scoring mode, either LOGARITHMIC (0), SCALE (1), ORIGINAL (2), RIBOSUM (3). (0)",
                                         ArgParseArgument::INTEGER, "INT"));

        addOption(parser, ArgParseOption("fsw", "fixedStructWeight",
                                         "define the weight of _half_ an interaction match for fixed structures",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("scal", "scalingFactor",
                                         "if structurescoring=SCALING then we have to give a scaling factor",
                                         ArgParseArgument::DOUBLE, "DOUBLE"));

        addOption(parser, ArgParseOption("m", "matching",
                                         "Lookahead for greedy matching algorithm. Value 0 uses LEMON instead. (5)",
                                         ArgParseArgument::INTEGER, "INT"));

        addOption(parser, ArgParseOption("j", "numThreads",
                                         "Use the number of specified threads. Decrease if memory problems occur. (5)",
                                         ArgParseArgument::INTEGER, "INT"));

        ArgumentParser::ParseResult parseResult = parse(parser, argc, argv);
        if (parseResult != ArgumentParser::PARSE_OK)
            return parseResult == ArgumentParser::ParseResult::PARSE_ERROR ? Status::EXIT_ERROR : Status::EXIT_OK;

        getOptionValue(_VERBOSE_LEVEL, parser, "verbose");
        getOptionValue(affineLinearDgs, parser, "affineLinearDgs");
        getOptionValue(alignLocally, parser, "local");
        getOptionValue(thrBppm, parser, "thrBppm");
        getOptionValue(numIterations, parser, "iterations");
        getOptionValue(maxNondecrIterations, parser, "maxNondecreasingIterations");
        getOptionValue(epsilon, parser, "epsilon");
        getOptionValue(stepSizeFactor, parser, "stepSizeFactor");
        getOptionValue(laraScoreMatrixName, parser, "laraScoreMatrixName");
        getOptionValue(suboptimalDiff, parser, "suboptimalDiff");
        getOptionValue(laraGapOpen, parser, "laraGapOpen");
        getOptionValue(laraGapExtend, parser, "laraGapExtend");
        getOptionValue(sequenceScale, parser, "sequenceScale");
        getOptionValue(balance, parser, "balance");
        getOptionValue(structureScoring, parser, "structureScoring");
        getOptionValue(fixedStructWeight, parser, "fixedStructWeight");
        getOptionValue(scalingFactor, parser, "scalingFactor");
        getOptionValue(tcoffeeLocation, parser, "tcoffeeLocation");
        getOptionValue(tcoffeeLibMode, parser, "tcoffeeLibMode");
        getOptionValue(inFileRef, parser, "inFileRef");
        getOptionValue(matching, parser, "matching");
        getOptionValue(num_threads, parser, "numThreads");

        getOptionValue(inFile, parser, "inFile");
        unsigned numDotplots = getOptionValueCount(parser, "dotplots");
        dotplotFile.resize(numDotplots);
        for (unsigned idx = 0; idx < numDotplots; ++idx)
        {
            getOptionValue(dotplotFile[idx], parser, "dotplots", idx);
        }

        if (empty(inFile) && dotplotFile.empty())
        {
//            std::cerr << "ERROR: You have to specify either -i or -d to pass input files." << std::endl;
            printShortHelp(parser);
            return Status::EXIT_ERROR;
        }

        if (num_threads == 0u)
        {
            unsigned nthreads = std::thread::hardware_concurrency() / 2u;
            if (nthreads != 0)
                num_threads = nthreads;
            else
                num_threads = 1u;
        }

        getOptionValue(outFile, parser, "outFile");
        if (isSet(parser, "outFile"))
        {
            _LOG(1, "The specified output file is " << outFile << std::endl);
        }
        else
        {
            CharString tmpDirectory;
            getOptionValue(tmpDirectory, parser, "tmpDir");
            if (!isSet(parser, "tmpDir"))
            {
                tmpDirectory = SEQAN_TEMP_FILENAME();
                // remove "/test_file" suffix
                erase(tmpDirectory, length(tmpDirectory) - 10u, length(tmpDirectory));
            }
            tmpDir = tmpDirectory;
            _LOG(2, "The absolute path where to create the tmpDir is " << tmpDir << std::endl);
        }

        // set score matrix
        laraScoreMatrix.data_gap_extend = laraGapExtend;
        laraScoreMatrix.data_gap_open   = laraGapOpen;
        if (empty(laraScoreMatrixName))
        {
            _LOG(2, "Predefined Ribosum65 matrix will be used." << std::endl);
            setRnaScoreMatrix(laraScoreMatrix, Ribosum65N());
        }
        else if (loadScoreMatrix(laraScoreMatrix, toCString(laraScoreMatrixName)))
        {
            _LOG(2, "Provided scoring matrix will be used: " << laraScoreMatrixName << std::endl);
        }
        else
        {
            std::cerr << "Matrix file could not be opened: " << laraScoreMatrixName
                      << "). Predefined Ribosum65 matrix will be used." << std::endl;
            setRnaScoreMatrix(laraScoreMatrix, Ribosum65N());
        }

        return Status::CONTINUE;
    }
};

} // namespace lara
