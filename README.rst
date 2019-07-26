LaRA 2: Lagrangian Relaxed structural Alignment
==============================================

LaRA 2 is an improved version of LaRA, a tool for sequence-structure alignment of RNA sequences. It...

* employs methods from combinatorial optimization to compute feasible solutions for an integer linear program
* computes all pairwise sequence-structure alignments of the input sequences and passes this information on to
  T-Coffee which computes a multiple sequence-structure alignment given the pairwise alignments
* can read many input formats for RNA structure, e.g. Dot-bracket notation, Stockholm, Vienna format
* is implemented to use multiple threads on your machine and runs therefore very fast
* is based on the SeqAn library, currently version 2
* is well-documented and easy to use


Download instructions
---------------------

Clone the repository and use the `--recursive` option for downloading SeqAn 2.4 and Lemon 1.3.1 as submodules.

::

% git clone --recursive https://github.com/seqan/lara.git

Alternatively, you can download a zip package of the repository via the green button at the top of the github page.
If you do so, please unzip the file into a new subdirectory `lara` and download the dependencies separately.


Requirements
------------

* platforms: Linux, MacOS
* compiler: gcc ≥ 5 or clang ≥ 3.8 or icc ≥ 17
* cmake ≥ 3

LaRA is dependent on the following libraries:

* `SeqAn 2.4 <https://github.com/seqan/seqan.git>`__
* `Lemon 1.3.1 <https://github.com/seqan/lemon.git>`__

Optionally, LaRA can predict the RNA structures for you if you provide

* `ViennaRNA 2 <https://www.tbi.univie.ac.at/RNA/>`__


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

The pairwise structural alignments are printed to stdout in the T-Coffee Library format.
If you want to store the result in a file, please use the `-w` option or redirect the output.

::

% bin/lara -i sequences.fasta -w results.lib
% bin/lara -i sequences.fasta  > results.lib

We recommend you to specify the number of threads with the `-j` option, e.g. to execute 4 alignments in parallel.
If you specify `-j 0` the program tries to detect the maximal number of threads available on your machine.

::

% bin/lara -i sequences.fasta -j 4

For a list of options, please see the help message:

::

% bin/lara --help


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
