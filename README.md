MXfold: the max-margin based RNA folding algorithm
=========================================================================

Requirements
---------------

* C++11 compatible compiler (tested on Apple LLVM version 6.1.0 and GCC version 4.8.1)
* [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/) (>= 2.3)

Install
---------

	./configure --with-vienna-rna=/path/to/vienna-rna
    make

or

    export PKG_CONFIG_PATH=/path/to/vienna-rna/lib/pkgconfig:${PKG_CONFIG_PATH}
    cmake -DCMAKE_BUILD_TYPE=Release .
    make


Usage
------

MXfold can take a FASTA formatted RNA sequence as input, then predicts
its secondary structure.

	% mxfold test.fa
	> DS4440
    GGAUGGAUGUCUGAGCGGUUGAAAGAGUCGGUCUUGAAAACCGAAGUAUUGAUAGGAAUACCGGGGGUUCGAAUCCCUCUCCAUCCG
    >structure
    (((((((........(((((..(((.......)))...)))))..(((((......))))).(((((.......)))))))))))).

License
---------

Copyright (c) 2017 Kengo Sato, Manato Akiyama  
Released under the MIT license  
http://opensource.org/licenses/mit-license.php


Acknowledgments
--------------------

MXfold is based on the source code of [CONTRAfold](http://contra.stanford.edu/contrafold/).


References
-------------

* Akiyama, M., Sato, K., Sakakibara, Y.: A max-margin training of RNA
  secondary structure prediction integrated with the thermodynamic
  model, submitted. [preprint](https://www.biorxiv.org/content/early/2017/10/18/205047)
