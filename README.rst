LaRA2: Lagrangian Relaxed structural Alignment
-----------------------------------------------------

LaRA2 is an improved version of LaRA, a tool for sequence-structure alignment of RNA sequences. It...

* employs methods from combinatorial optimization to compute feasible solutions for an integer linear program
* computes all pairwise sequence-structure alignments of the input sequences and passes this information on to T-Coffee which computes a multiple sequence-structure alignment given the pairwise alignments
* can read many input formats for RNA structure, e.g. Dot-bracket notation, Stockholm, Vienna format
* is based on SeqAn2
* is well-documented and easy to use

downloads and installation
--------------------------

+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
|  **Executables**                                                                                                                                                                                                                      |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.disk.download.png        | will arrive here when the first release is ready                                                                   |
|    :alt: Download Executables                                                                                    |                                                                                                                    |
|    :target: https://github.com/seqan/lara/releases                                                               |                                                                                                                    |
|    :width: 76px                                                                                                  |                                                                                                                    |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
|  **Source code**                                                                                                                                                                                                                      |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.column.three.png         | You can build LaRA2 from source which will result in binaries optimized for your                                   |
|    :alt: Build from source                                                                                       | specific system (and thus faster). For instructions, please see below.                                             |
|    :target: https://github.com/seqan/lara/wiki                                                                   |                                                                                                                    |
|    :width: 76px                                                                                                  |                                                                                                                    |
+------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

download and usage instructions
------------------


Clone the repository:

::

% git clone --recursive https://github.com/seqan/lara.git

::

For a list of options, see the help page:

::

% bin/lara --help

authorship and copyright
------------------------

LaRA2 is being developed by `Jörg Winkler <mailto:j.winkler@fu-berlin.de>`__ and `Gianvito Urgese <mailto:gianvito.urgese@polito.it>`__, but it incorporates a lot of work from other members of the `SeqAn project <http://www.seqan.de>`__.

feedback & updates
------------------

+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.social.github.octocat.png | You can ask questions and report bugs on the `github tracker <https://github.com/seqan/lara/issues>`__ .           |
|    :alt: GitHub                                                                                                   | Please also `subscribe <https://github.com/seqan/lara/subscription>`__ and/or star us!                             |
|    :target: https://github.com/seqan/lara/issues                                                                  |                                                                                                                    |
|    :width: 76px                                                                                                   |                                                                                                                    |
+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+
| .. image:: https://raw.githubusercontent.com/seqan/lambda/gh-pages/images_readme/appbar.social.twitter.png        | You can also follow SeqAn on `twitter <https://twitter.com/SeqAnLib>`__ to receive updates on LaRA.                |
|    :alt: Newsletter                                                                                               |                                                                                                                    |
|    :target: https://twitter.com/SeqAnLib                                                                          |                                                                                                                    |
|    :width: 76px                                                                                                   |                                                                                                                    |
+-------------------------------------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------------+

*icons on this page by Austin Andrews / https://github.com/Templarian/WindowsIcons*
