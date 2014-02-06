TargetFinder
============

Plant small RNA target prediction tool

About
============
TargetFinder will computationally predict small RNA binding sites on target transcripts from a sequence database. This is done by aligning the input small RNA sequence against all transcripts, followed by site scoring using a position-weighted scoring matrix.

Installing TargetFinder
============

Minimum Requirements
------------

* Required Software Packages:
  * Perl (v5.8+)
  * FASTA34 binaries from (http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml)

Installation
------------

1.  Make sure that 'mkdir' is in your path.
2.  Set an environmental variable 'HOME' as the location of your home directory. Note: you must have read/write access to this directory so that a temporary directory can be made.
3.  Set an environmental variable 'FASTA34' as the location of the fasta34 binary.
4.  Make sure that the required modules are in your Perl library path.
5.  Make sure that targetfinder.pl is executable (use 'chmod' to change privileges).

Using TargetFinder
============

Required arguments
------------
  -s     Small RNA sequence (either RNA or DNA).
  
  -d     FASTA-formated sequence file containing potential target sequences.

Optional Arguments
------------
  -q     Query sequence name (DEFAULT = query).
  
  -c     Prediction score cutoff value (DEFAULT = 4).
  
  -h     Shows the help menu.

Output
------------
targetfinder.pl writes all output to the terminal (STOUT). To save the output to a file use '>' to redirect output to a file.

  ex.

  ./targetfinder.pl -s UGUGUUCUCAGGUCACCCCUU -d arab_cdna -q miR399a > miR399a_predicted_targets.txt

Each predicted target site is printed out separately. The output consists of two parts. The first is a description line and the second is a base-pairing diagram of the target and small RNA (query) sequence. The description line contains the query name (query=name), the description line from the target sequence database (target=target description), and the target prediction score (score=prediction score).

  ex.
 
  query=miR399a, target=AT2G33770.1 | Symbol: None |  ubiquitin-conjugating enzyme family protein, low similarity to u, score=1.5

The base-pairing diagram has the target site sequence on top in 5'-3' orientation and the query sequence on the bottom in 3'-5' orientation. Between the target site sequece and the query sequence are base pair
symbols. A ":" (colon) symbol represents an ordinary Watson-Crick base pair, a "." (period) represents a G:U base pair, and a " " (space) represents a mismatch.

  ex.

` 
target  5' UAGGGCAAAUCUUCUUUGGCA 3'
           .:::::::::::.::::::::
query   3' GUCCCGUUUAGAGGAAACCGU 5'
`  
 
If a small RNA is predicted to target a sequence more than once, each target site will be output as separate output. Below is an example of output for miR399a and its target AT2G33770. miR399a has five
target sites in the 5'UTR of AT2G33770.

  query=miR399a, target=AT2G33770.1 | Symbol: None |  ubiquitin-conjugating enzyme family protein, low similarity to u, score=1.5
  
  target  5' UAGGGCAAAUCUUCUUUGGCA 3'
             .:::::::::::.::::::::
  query   3' GUCCCGUUUAGAGGAAACCGU 5'
  
  query=miR399a, target=AT2G33770.1 | Symbol: None |  ubiquitin-conjugating enzyme family protein, low similarity to u, score=1.5
  
  target  5' UAGGGCAUAUCUCCUUUGGCA 3'
             .:::::: :::::::::::::
  query   3' GUCCCGUUUAGAGGAAACCGU 5'
  
  query=miR399a, target=AT2G33770.1 | Symbol: None |  ubiquitin-conjugating enzyme family protein, low similarity to u, score=1.5
  
  target  5' UAGAGCAAAUCUCCUUUGGCA 3'
             .:: :::::::::::::::::
  query   3' GUCCCGUUUAGAGGAAACCGU 5'
  
  query=miR399a, target=AT2G33770.1 | Symbol: None |  ubiquitin-conjugating enzyme family protein, low similarity to u, score=1.5
  
  target  5' UUGGGCAAAUCUCCUUUGGCA 3'
             . :::::::::::::::::::
  query   3' GUCCCGUUUAGAGGAAACCGU 5'
  
  query=miR399a, target=AT2G33770.1 | Symbol: None |  ubiquitin-conjugating enzyme family protein, low similarity to u, score=2.5
  
  target  5' UCGAGCAAAUCUCCUUUGGCA 3'
             . : :::::::::::::::::
  query   3' GUCCCGUUUAGAGGAAACCGU 5'

Method
============
targetfinder.pl searches for potential miRNA target sites in a FASTA-formated sequence database using three main steps.

1.  The small RNA query sequence is aligned to every sequence in the FASTA-formated sequence database using the alignment program FASTA34.
2.  The FASTA34 alignments are converted into RNA duplexes.
3.  Each duplex is scored using a position-dependent scoring matrix.

FASTA34 is used to identify the best complementary regions between the small RNA query sequence and every sequence in the FASTA-formated sequence database. This script runs FASTA34 with the following settings:

  -n     Forces the small RNA query sequence to be treated as nucleotide sequence.
  -H     Suppresses the normal histogram output of FASTA34.
  -Q     Runs FASTA34 in "quiet" mode.
  -f     Gap opening penalty (set to -16).
  -g     Gap extention penalty (set to -10).
  -r     Match reward/mismatch penalty (set to +15/-10).
  -w     Alignment output line length (set to 100).
  -W     Additional sequence context in the output (set to 25).
  -b     The number of alignments to include in the output (set to 100000).
  -i     Limits FASTA34 alignments to reverse complement matches only.
  -U     Changes scoring matrix to allow for G:A, T:C, or U:C matches.
  ktup   Word size for seed matches that FASTA34 uses to build alignments (set to 1).

FASTA34 output is read directly into this script.  Each alignment is converted to a RNA duplex by complementing the small RNA query sequence. Each RNA duplex is scored using the following scoring metric and rule set:

1.  Mismatches, single-nucleotide gaps or single-nucleotide bulges are assesed a penalty of +1.
2.  G:U base pairs are assessed a penalty of +0.5.
3.  Penalty scores are doubled at positions 2-13 relative to the 5' end of the small RNA query sequence.
4.  Duplexes are rejected if they:
  * have more than one single-nucleotide bulge or gap.
  * have more than seven total mismatches, G:U base pairs, bulges and gaps.
  * have more than four total mismatches or four total G:U base pairs.

Predicted targets are printed out if they are equal to or lower than the cutoff score specified.

References
============
1.  Allen E, Xie Z, Gustafson AM, Carrington JC (2005) microRNA-directed phasing during trans-acting siRNA biogenesis in plants. Cell 121: 207Ð221. doi:[10.1016/j.cell.2005.04.004](http://dx.doi.org/10.1016/j.cell.2005.04.004).
2.  Fahlgren N, Howell MD, Kasschau KD, Chapman EJ, Sullivan CM, Cumbie JS, Givan SA, Law TF, Grant SR, Dangl JL, Carrington JC (2007) High-throughput sequencing of *Arabidopsis* microRNAs: Evidence for frequent birth and death of *MIRNA* genes. PLoS ONE 2: e219. doi:[10.1371/journal.pone.0000219](http://dx.doi.org/10.1371/journal.pone.0000219).