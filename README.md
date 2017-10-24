MXfold: the max-margin based RNA folding algorithm
=========================================================================

Requirements
------------

* C++11 compatible compiler (tested on Apple LLVM version 6.1.0 and GCC version 4.8.1)
* [Vienna RNA package](https://www.tbi.univie.ac.at/RNA/) (>= 2.3)

Install
-------

	./configure --with-vienna-rna=/path/to/vienna-rna
    make
    make install # optional

Usage
-----

MXfold can take a FASTA formatted RNA sequence as input, then predicts
its secondary structure.

	% mxfold test.fa
	> DS4440
    GGAUGGAUGUCUGAGCGGUUGAAAGAGUCGGUCUUGAAAACCGAAGUAUUGAUAGGAAUACCGGGGGUUCGAAUCCCUCUCCAUCCG
    >structure
    (((((((........(((((..(((.......)))...)))))..(((((......))))).(((((.......)))))))))))).

References
----------

* Akiyama, M., Sato, K., Sakakibara, Y.: A max-margin training of RNA
  secondary structure prediction integrated with the thermodynamic
  model, submitted. [preprint](https://www.biorxiv.org/content/early/2017/10/18/205047)
