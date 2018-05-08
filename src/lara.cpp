#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/rna_io.h>

extern "C" {
    #include <ViennaRNA/data_structures.h>
    #include <ViennaRNA/params.h>
    #include <ViennaRNA/utils.h>
    #include <ViennaRNA/eval.h>
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/part_func.h>
    #include <ViennaRNA/PS_dot.h>
}

using namespace seqan;

int main (int argc, char const ** argv)
{
    std::cout << "This is an initial set-up for Lara2. See development branch on github for current implementation "
                 "status." << std::endl << argv[0];
    for (int i = 1; i <= argc; ++i)
    {
        std::cout << " " << argv[i];
    }
    std::cout << std::endl;
}
