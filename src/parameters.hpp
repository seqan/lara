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

/*!\file parameters.hpp
 * \brief This file contains all settings and parameters of LaRA.
 */

#include <seqan/arg_parse.h>

#include "data_types.hpp"

namespace lara
{

struct Parameters
{
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
    double                   thrBppm{1e-15}; // old Lara: 0.1
    // number of iterations
    unsigned                 iterations{500u};
    // number of non-decreasing iterations
    unsigned                 nonDecreasingIterations{50u};
    // value to be considered for the equality of upper and lower bounds difference
    double                   epsilon{0.0001};
    // my, necessary for computing appropriate step sizes
    double                   stepSizeScaling{1.0};
    // scoring matrix name that should be used for scoring alignment edges in the actual problem
    seqan::CharString        laraScoreMatrixName{};
    //    Score<double, ScoreMatrix<Rna5, Default> > laraScoreMatrix;
    RnaScoreMatrix           laraScoreMatrix;
    //    TScoringSchemeRib laraScoreMatrixRib;
    // Gap open and extend costs for generating the alignment edges
    double                   generatorGapOpen{-6.0};
    double                   generatorGapExtend{-2.0};
    // Gap open and extend costs for generating the alignment edges
    double                   laraGapOpen{-6.0};
    double                   laraGapExtend{-2.0};
    // scaling factor for the scores of the alignment edges
    // Specifies the contribution of the sequence scores (specified by the larascore matrix) to the overall structural
    // alignment.
    double                   sequenceScale{1.0};
    // gap penalty for RSA
    //    double rsaGapPenalty{3.0};
    // scoring mode, either LOGARITHMIC, SCALE, ORIGINAL, RIBOSUM
    unsigned                 structureScoring{ScoringMode::LOGARITHMIC};
    // define the weight of _half_ an interaction match for fixed structures
    double                   fixedStructWeight{8.0};
    // if structureScoring=SCALING then we have to give a scaling factor
    double                   scalingFactor{1.0};
    // specify the location of T-COFFEE
    seqan::CharString        tcoffeeLocation{"t_coffee/t_coffee_5.05"};
    // specify the method to be used to create the T-Coffe library
    unsigned                 tcoffeLibMode{TCoffeeMode::SWITCH};
    // verbosity level (0-3)
    unsigned                 verbose{0u};
};

inline int setParameters(Parameters & params, int argc, char const **argv)
{
    using namespace seqan;
    ArgumentParser parser;

    setAppName(parser, "LaRA");
    setShortDescription(parser, "Lagrangian Relaxed structural Alignment");
    setCategory(parser, "RNA structural alignment algorithm");
    setVersion(parser, "2.0");
    setDate(parser, "2018");
    //setDateAndVersion(parser);
    //setDescription(parser);

    addUsageLine(parser, "./lara <\\fI-i inFile\\fP> [\\fI-w outFile\\fP] [\\fI -parameters\\fP]");

    addOption(parser, ArgParseOption("v", "verbose",
                                     "0: no additional outputs, 1: global statistics, "
                                     "2: extensive statistics for each batch of reads, 3: Debug output. (1)",
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
                                     "Use dotplot files (.ps) as bpp input",
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

    addOption(parser, ArgParseOption("tcm", "tcoffeLibMode",
                                     "method used to create the T-Coffe library either 0: PROPORTIONAL, 1: SWITCH, "
                                     "2: ALLINTER, 3: FIXEDINTER. (0)",
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

    addOption(parser, ArgParseOption("nditer", "nonDecreasingIterations",
                                     "number of non-decreasing iterations. (50)",
                                     ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("ep", "epsilon",
                                     "value to be considered for the equality of upper and lower bounds difference",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("my", "stepSizeScaling",
                                     "necessary for computing appropriate step sizes.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("lsm", "laraScoreMatrixName",
                                     "scoring matrix name that should be used for scoring alignment edges in the "
                                     "actual problem",
                                     ArgParseOption::STRING));

    addOption(parser, ArgParseOption("ggo", "generatorGapOpen",
                                     "Gap open costs for generating the alignment edges.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("gge", "generatorGapExtend",
                                     "Gap extend costs for generating the alignment edges.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("lgo", "laraGapOpen",
                                     "Gap open costs for generating the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("lge", "laraGapExtend",
                                     "Gap extend costs for generating the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("ssc", "sequenceScale",
                                     "scaling factor for the scores of the alignment edges",
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

    ArgumentParser::ParseResult parseResult = parse(parser, argc, argv);
    if (parseResult != ArgumentParser::PARSE_OK)
        return parseResult == ArgumentParser::ParseResult::PARSE_ERROR ? 1 : 0;

    getOptionValue(params.verbose, parser, "verbose");
    getOptionValue(params.affineLinearDgs, parser, "affineLinearDgs");
    getOptionValue(params.alignLocally, parser, "local");
    getOptionValue(params.thrBppm, parser, "thrBppm");
    getOptionValue(params.iterations, parser, "iterations");
    getOptionValue(params.nonDecreasingIterations, parser, "nonDecreasingIterations");
    getOptionValue(params.epsilon, parser, "epsilon");
    getOptionValue(params.stepSizeScaling, parser, "stepSizeScaling");
    getOptionValue(params.laraScoreMatrixName, parser, "laraScoreMatrixName");
    getOptionValue(params.generatorGapOpen, parser, "generatorGapOpen");
    getOptionValue(params.generatorGapExtend, parser, "generatorGapExtend");
    getOptionValue(params.laraGapOpen, parser, "laraGapOpen");
    getOptionValue(params.laraGapExtend, parser, "laraGapExtend");
    getOptionValue(params.sequenceScale, parser, "sequenceScale");
    getOptionValue(params.structureScoring, parser, "structureScoring");
    getOptionValue(params.fixedStructWeight, parser, "fixedStructWeight");
    getOptionValue(params.scalingFactor, parser, "scalingFactor");
    getOptionValue(params.tcoffeeLocation, parser, "tcoffeeLocation");
    getOptionValue(params.tcoffeLibMode, parser, "tcoffeLibMode");
    getOptionValue(params.inFileRef, parser, "inFileRef");

    getOptionValue(params.inFile, parser, "inFile");
    if (empty(params.inFile))
        return 1;

    unsigned numDotplots = getOptionValueCount(parser, "dotplots");
    params.dotplotFile.resize(numDotplots);
    for (unsigned idx = 0; idx < numDotplots; ++idx)
    {
        getOptionValue(params.dotplotFile[idx], parser, "dotplots", idx);
    }

    getOptionValue(params.outFile, parser, "outFile");
    if (isSet(parser, "outFile"))
    {
        _V(params, "The specified output file is " << params.outFile);
    }
    else
    {
        CharString tmpDir;
        getOptionValue(tmpDir, parser, "tmpDir");
        if (!isSet(parser, "tmpDir"))
        {
            tmpDir = SEQAN_TEMP_FILENAME();
            // remove "/test_file" suffix
            erase(tmpDir, length(tmpDir) - 10u, length(tmpDir));
        }
        params.tmpDir = tmpDir;
        _V(params, "The absolute path where to create the tmpDir is " << tmpDir);
    }

    return 2;
}

} // namespace lara

