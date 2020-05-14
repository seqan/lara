What LaRA 2 does for you
========================

LaRA 2 is an improved version of LaRA, a tool for sequence-structure alignment of RNA sequences. It...

* computes all pairwise sequence-structure alignments of the input sequences
* produces files that can be processed with T-Coffee or MAFFT to compute a multiple sequence-structure alignment
* employs methods from combinatorial optimization to compute feasible solutions for an integer linear program
* can read many input formats for RNA structure, e.g. Dot-bracket notation, Stockholm, Vienna format
* is implemented to use multiple threads on your machine and runs therefore very fast
* has a vectorized alignment kernel, which computes the results even faster
* is based on the SeqAn library, currently version 2
* is well-documented and easy to use


Download instructions
---------------------

Clone the repository and use the *-\-recurse-submodules* option for downloading SeqAn and Lemon as submodules.

```commandline
git clone --recurse-submodules https://github.com/seqan/lara.git
```

Alternatively, you can download a package of the repository via the buttons at the top of this page.
If you do so, please extract the file into a new subdirectory named lara and download the dependencies separately.


Requirements
------------

* platforms: Linux, MacOS
* compiler: gcc ≥ 5 or clang ≥ 3.8 or icc ≥ 17
* cmake ≥ 3.8

LaRA is dependent on the following libraries:

* [SeqAn 2.4](https://github.com/seqan/seqan.git)
* [Lemon 1.3.1](https://github.com/seqan/lemon.git)

If you have not performed a recursive clone above, simply run the following command in the lara directory
to download them.

```commandline
cd lara
git submodule update --init --recursive
cd ..
```

Optionally, LaRA can predict the RNA structures for you if you provide

* [ViennaRNA 2](https://www.tbi.univie.ac.at/RNA/)

To process the output for multiple alignments (3 or more sequences), you need either

* [T-Coffee 13](https://github.com/cbcrg/tcoffee) or
* [MAFFT 7.453 for LaRA](https://github.com/bioinformatics-polito/LaRA2-mafft)


Build instructions
------------------

Please create a new directory and build the program for your platform.

```commandline
mkdir build
cd build
cmake ../lara     # specify the path to the lara directory
make
cd ..
```


First steps to use LaRA 2
=========================

After building the program binary, running LaRA is as simple as

```commandline
build/lara -i sequences.fasta
```

Note that for passing sequence files you need the ViennaRNA dependency, as the program must predict structures.
Instead, you can pass at least two dot plot files, which contain the base pair probabilities for a single sequence each.

```commandline
build/lara -d seq1_dp.ps -d seq2_dp.ps
```

The pairwise structural alignments are printed to stdout in the T-Coffee Library format (see below).
If you want to store the result in a file, please use the *-w* option or redirect the output.

```commandline
build/lara -i sequences.fasta -w results.lib
build/lara -i sequences.fasta  > results.lib
```

We recommend you to specify the number of threads with the *-j* option, to execute for instance 4 alignments in
parallel. If you specify *-j 0* the program tries to detect the maximal number of threads available on your machine.

```commandline
build/lara -i sequences.fasta -j 4
```

For a list of options, please see the help message:

```commandline
build/lara --help
```


Output format
-------------

The result of LaRA is a T-Coffee library file and its format is documented
[here](http://www.tcoffee.org/Projects/tcoffee/documentation/index.html#t-coffee-lib-format-01).
It contains the structural scores for each residue pair of each computed sequence pair.
This file is the input for T-Coffee, which computes the multiple alignment based on the scores:

```commandline
t_coffee -lib results.lib
```


Authorship & Copyright
----------------------

LaRA 2 is being developed by [Jörg Winkler](mailto:j.winkler@fu-berlin.de) and
[Gianvito Urgese](mailto:gianvito.urgese@polito.it), but it incorporates a lot of work
from other members of the [SeqAn project](http://www.seqan.de).

