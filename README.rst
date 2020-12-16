LaRA 2: Lagrangian Relaxed structural Alignment
===============================================

LaRA 2 is an improved version of LaRA, a tool for sequence-structure alignment of RNA sequences. It...

* computes all pairwise sequence-structure alignments of the input sequences
* produces files that can be processed with T-Coffee or MAFFT to compute a multiple sequence-structure alignment
* employs methods from combinatorial optimization to compute feasible solutions for an integer linear program
* can read many input formats for RNA structure, e.g. Dot-bracket notation, Stockholm, Vienna format
* is implemented to use multiple threads on your machine and runs therefore very fast
* has a vectorized alignment kernel, which computes the results even faster
* is based on the SeqAn library, currently version 2
* is well-`documented <https://seqan.github.io/lara/>`__ and easy to use


Download instructions
---------------------

Clone the repository and use the *-\-recurse-submodules* option for downloading SeqAn and Lemon as submodules.

::

  % git clone --recurse-submodules https://github.com/seqan/lara.git

Alternatively, you can download a zip package of the repository via the green button at the top of the github page.
If you do so, please unzip the file into a new subdirectory named *lara* and download the dependencies separately.


Requirements
------------

* platforms: Linux, MacOS
* compiler: gcc ≥ 5 or clang ≥ 3.8 or icc ≥ 17
* cmake ≥ 3.8

LaRA is dependent on the following libraries:

* `SeqAn 2.4 <https://github.com/seqan/seqan.git>`__
* `Lemon 1.3.1 <https://github.com/seqan/lemon.git>`__

To process the output for multiple alignments (3 or more sequences), you need either

* `T-Coffee 13 <https://github.com/cbcrg/tcoffee>`__ or
* `MAFFT 7.453 for LaRA <https://github.com/bioinformatics-polito/LaRA2-mafft>`__

Optionally, LaRA can predict the RNA structures for you if you provide

* `ViennaRNA 2 <https://www.tbi.univie.ac.at/RNA/>`__

*Note:* Users reported problems with installing ViennaRNA, so we provide some hints here.

1. Install the `GNU MPFR Library <https://www.mpfr.org/>`__ first.
2. Exclude unnecessary components of ViennaRNA:
   ``./configure --without-swig --without-kinfold --without-forester --without-rnalocmin --without-gsl``
3. If you have linker issues use
   ``./configure --disable-lto``
4. If your system supports SSE4.1 instructions then we recommend
   ``./configure --enable-sse``

If you have further suggestions, we are happy to add them here.


Build instructions
------------------

Please create a new directory and build the program for your platform.

::

  % mkdir bin
  % cd bin
  % cmake ../lara
  % make
  % cd ..


Usage
-----

After building the program binary, running LaRA is as simple as

::

  % bin/lara -i sequences.fasta

Note that for passing sequence files you need the ViennaRNA dependency, as the program must predict structures.
Instead, you can pass at least two dot plot files, which contain the base pair probabilities for a single sequence each.

::

  % bin/lara -d seq1_dp.ps -d seq2_dp.ps

The pairwise structural alignments are printed to stdout in the T-Coffee Library format (see below).
If you want to store the result in a file, please use the *-w* option or redirect the output.

::

  % bin/lara -i sequences.fasta -w results.lib
  % bin/lara -i sequences.fasta  > results.lib

We recommend you to specify the number of threads with the *-j* option, e.g. to execute 4 alignments in parallel.
If you specify *-j 0* the program tries to detect the maximal number of threads available on your machine.

::

  % bin/lara -i sequences.fasta -j 4

For a list of options, please see the help message:

::

  % bin/lara --help


Output format
-------------

Each output format is sorted primarily by the first and subsequently by the second sequence index.

for multiple alignments with T-Coffee
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The result of LaRA is a T-Coffee library file and its format is documented
`here <http://www.tcoffee.org/Projects/tcoffee/documentation/index.html#t-coffee-lib-format-01>`__.
It contains the structural scores for each residue pair of each computed sequence pair.
This file is the input for T-Coffee, which computes the multiple alignment based on the scores:

::

  % bin/t_coffee -lib results.lib

for multiple alignments with MAFFT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LaRA has an additional output format that can be read by the MAFFT framework.
Each pairwise alignment produces three lines:
a description line composed of the two sequence ids and the two gapped sequences of the alignment.

::

  > first id && second id
  AACCG-UU
  -ACCGGUU
  > first id && third id
  AA-CCGUU
  AAGCCGUU

MAFFT invokes LaRA with the option *-o pairs* for receiving this output format.

for pairwise alignments
~~~~~~~~~~~~~~~~~~~~~~~

LaRA can produce the aligned FastA format, which is recommended for a single pairwise alignment.
It looks like a normal FastA file with gap symbols in the sequences:

::

  > first id
  AACCG-UU
  > second id
  -ACCGGUU

You need to pass the option *-o fasta* to the LaRA call for getting this output format.

LaRA prints a warning if you use this format with more than two sequences.
Using this format with 3 or more sequences is possible but not recommended, because additional pairwise alignments
will simply be appended to the file, and it may be hard to distinguish the pairs later.
In addition, this can confuse other programs, which expect a single multiple sequence alignment
as produced by MAFFT or T-Coffee.


Authorship & Copyright
----------------------

LaRA 2 is being developed by `Jörg Winkler <mailto:j.winkler@fu-berlin.de>`__ and
`Gianvito Urgese <mailto:gianvito.urgese@polito.it>`__, but it incorporates a lot of work
from other members of the `SeqAn project <http://www.seqan.de>`__.


Feedback & Updates
------------------

+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.social.github.octocat.png | You can ask questions and report bugs on the `github tracker <https://github.com/seqan/lara/issues>`__.            |
|    :alt: GitHub                                                                                                   | Please also `subscribe <https://github.com/seqan/lara/subscription>`__ and/or star us!                             |
|    :target: https://github.com/seqan/lara/issues                                                                  |                                                                                                                    |
|    :width: 76px                                                                                                   |                                                                                                                    |
+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.social.twitter.png        | You can also follow SeqAn on `twitter <https://twitter.com/SeqAnLib>`__ to receive updates on LaRA.                |
|    :alt: Newsletter                                                                                               |                                                                                                                    |
|    :target: https://twitter.com/SeqAnLib                                                                          |                                                                                                                    |
|    :width: 76px                                                                                                   |                                                                                                                    |
+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

*Icons on this page by Austin Andrews: https://github.com/Templarian/WindowsIcons*
