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
  * FASTA35 binaries from (http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml)

Installation
------------

1.  TargetFinder uses the environmental variable TMPDIR as a temporary directory. If TMPDIR is undefined TargetFinder tries to use /tmp instead. If you want TargetFinder to write temporary files to a specific place, change the value of TMPDIR.
2.  Make sure that 'ssearch35_t' is in your path.
3.  Make sure that the required modules are in your Perl library path.
4.  Make sure that targetfinder.pl is executable (use 'chmod' to change privileges).

Using TargetFinder
============

Required arguments
------------
  -s <string>    Small RNA sequence (RNA or DNA, 5'->3').
  
  -d <file>      Target sequence database file (FASTA-format).

Optional Arguments
------------
  -q <string>    Query sequence name (DEFAULT = query).
  
  -c <float>     Prediction score cutoff value (DEFAULT = 4).
  
  -t <int>       Threads for parallel Smith-Waterman searches (DEFAULT = 1).

  -p <string>    Output format for small RNA-target pairs (DEFAULT = 'classic')
                 Available options: 'classic' (Original TargetFinder base-pairing format)
                                    'gff'     (Generic Feature Format)
                                    'json'    (JavaScript Object Notation)
                                    'table'   (Tab-deliminated Format)

  -r             Search reverse strand for targets?. Use this option if the database is genomic DNA.
  
  -h             Shows the help menu.

Output
------------
targetfinder.pl writes all output to the terminal (STOUT). To save the output to a file use '>' to redirect output to a file.

  example:

```
./targetfinder.pl -s UGCCAAAGGAGAUUUGCCCUG -d arab_cdna -q miR399a > miR399a_predicted_targets.txt
```

Each predicted target site is printed out separately. The output consists of two parts. The first is a description line and the second is a base-pairing diagram of the target and small RNA (query) sequence. The description line contains the query name (query=name), the description line from the target sequence database (target=target description), and the target prediction score (score=prediction score).

  example:

``` 
query=miR399a, target=AT2G33770.1 | Symbol: None |  ubiquitin-conjugating enzyme family protein, low similarity to u, score=1.5
```

The base-pairing diagram has the target site sequence on top in 5'-3' orientation and the query sequence on the bottom in 3'-5' orientation. Between the target site sequece and the query sequence are base pair
symbols. A ":" (colon) symbol represents an ordinary Watson-Crick base pair, a "." (period) represents a G:U base pair, and a " " (space) represents a mismatch.

  example:

```
target  5' UAGGGCAAAUCUUCUUUGGCA 3'  
           .:::::::::::.::::::::  
query   3' GUCCCGUUUAGAGGAAACCGU 5'
```
 
If a small RNA is predicted to target a sequence more than once, each target site will be output as separate output. Below is an example of output for miR399a and its target AT2G33770. miR399a has five
target sites in the 5'UTR of AT2G33770.

```
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
```

New Output Options
------------
In addition to the output described above ('classic' output), three new output format options were added to TargetFinder.

Generic Feature Format (GFF3):
```
./targetfinder.pl -s UGCCAAAGGAGAUUUGCCCUG -d arab_cdna -q miR399a -p gff > miR399a_predicted_targets.gff3

AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN  targetfinder  rna_target  607 627 1.5 + . smallRNA=miR399a;target_seq=UAGGGCAAAUCUUCUUUGGCA;base_pairs=.:::::::::::.::::::::;miR_seq=GUCCCGUUUAGAGGAAACCGU
AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN  targetfinder  rna_target  740 760 1.5 + . smallRNA=miR399a;target_seq=UAGGGCAUAUCUCCUUUGGCA;base_pairs=.:::::: :::::::::::::;miR_seq=GUCCCGUUUAGAGGAAACCGU
AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN  targetfinder  rna_target  829 849 1.5 + . smallRNA=miR399a;target_seq=UUGGGCAAAUCUCCUUUGGCA;base_pairs=. :::::::::::::::::::;miR_seq=GUCCCGUUUAGAGGAAACCGU
AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN  targetfinder  rna_target  943 963 1.5 + . smallRNA=miR399a;target_seq=UAGAGCAAAUCUCCUUUGGCA;base_pairs=.:: :::::::::::::::::;miR_seq=GUCCCGUUUAGAGGAAACCGU
AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN  targetfinder  rna_target  886 906 2.5 + . smallRNA=miR399a;target_seq=UCGAGCAAAUCUCCUUUGGCA;base_pairs=. : :::::::::::::::::;miR_seq=GUCCCGUUUAGAGGAAACCGU
```

Tab-deliminated Format:
```
miR399a	AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN	607	627	+	1.5	UAGGGCAAAUCUUCUUUGGCA	.:::::::::::.::::::::	GUCCCGUUUAGAGGAAACCGU
miR399a	AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN	740	760	+	1.5	UAGGGCAUAUCUCCUUUGGCA	.:::::: :::::::::::::	GUCCCGUUUAGAGGAAACCGU
miR399a	AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN	829	849	+	1.5	UUGGGCAAAUCUCCUUUGGCA	. :::::::::::::::::::	GUCCCGUUUAGAGGAAACCGU
miR399a	AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN	943	963	+	1.5	UAGAGCAAAUCUCCUUUGGCA	.:: :::::::::::::::::	GUCCCGUUUAGAGGAAACCGU
miR399a	AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN	886	906	+	2.5	UCGAGCAAAUCUCCUUUGGCA	. : :::::::::::::::::	GUCCCGUUUAGAGGAAACCGU
```

JavaScript Object Notation Format (JSON):
```
{
  "miR399a": {
    "hits" : [
      {
        "Target accession": "AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN",
        "Score": "1.5",
        "Coordinates": "607-627",
        "Strand": "+",
        "Target sequence": "UAGGGCAAAUCUUCUUUGGCA",
        "Base pairing": ".:::::::::::.::::::::",
        "amiRNA sequence": "GUCCCGUUUAGAGGAAACCGU"
      },
      {
        "Target accession": "AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN",
        "Score": "1.5",
        "Coordinates": "740-760",
        "Strand": "+",
        "Target sequence": "UAGGGCAUAUCUCCUUUGGCA",
        "Base pairing": ".::::::&nbsp:::::::::::::",
        "amiRNA sequence": "GUCCCGUUUAGAGGAAACCGU"
      },
      {
        "Target accession": "AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN",
        "Score": "1.5",
        "Coordinates": "829-849",
        "Strand": "+",
        "Target sequence": "UUGGGCAAAUCUCCUUUGGCA",
        "Base pairing": ".&nbsp:::::::::::::::::::",
        "amiRNA sequence": "GUCCCGUUUAGAGGAAACCGU"
      },
      {
        "Target accession": "AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN",
        "Score": "1.5",
        "Coordinates": "943-963",
        "Strand": "+",
        "Target sequence": "UAGAGCAAAUCUCCUUUGGCA",
        "Base pairing": ".::&nbsp:::::::::::::::::",
        "amiRNA sequence": "GUCCCGUUUAGAGGAAACCGU"
      },
      {
        "Target accession": "AT2G33770.1 | Symbols: UBC24, ATUBC24, PHO2 | phosphate 2 | chr2:14277558-14283040 REVERSE LEN",
        "Score": "2.5",
        "Coordinates": "886-906",
        "Strand": "+",
        "Target sequence": "UCGAGCAAAUCUCCUUUGGCA",
        "Base pairing": ".&nbsp:&nbsp:::::::::::::::::",
        "amiRNA sequence": "GUCCCGUUUAGAGGAAACCGU"
      }
    ]
  }
}
```

Method
============
targetfinder.pl searches for potential miRNA target sites in a FASTA-formated sequence database using three main steps.

1.  The small RNA query sequence is aligned to every sequence in the FASTA-formated sequence database using Smith-Waterman (SW) alignments implemented in the FASTA package (ssearch35_t).
2.  The SW alignments are converted into RNA duplexes.
3.  Each duplex is scored using a position-dependent scoring matrix.

SW alignments are used to identify the best complementary regions between the small RNA query sequence and every sequence in the FASTA-formated sequence database. This script runs ssearch35_t with the following settings:
  
  -n     Forces the small RNA query sequence to be treated as nucleotide sequence.  
  -H     Suppresses the normal histogram output.  
  -Q     Runs Smith-Waterman search in "quiet" mode.  
  -f     Gap opening penalty (set to -16).  
  -g     Gap extention penalty (set to -10).  
  -r     Match reward/mismatch penalty (set to +15/-10).  
  -w     Alignment output line length (set to 100).  
  -W     Additional sequence context in the output (set to 25).  
  -E     The E-value cutoff (set to 100000).  
  -i     Limits SW alignments to reverse complement matches only.  
  -U     Changes scoring matrix to allow for G:A, T:C, or U:C matches.  
  ktup   Word size for seed matches used to build alignments (set to 1).  

SW output is read directly into this script.  Each alignment is converted to a RNA duplex by complementing the small RNA query sequence. Each RNA duplex is scored using the following scoring metric and rule set:

1.  Mismatches, single-nucleotide gaps or single-nucleotide bulges are assesed a penalty of +1.
2.  G:U base pairs are assessed a penalty of +0.5.
3.  Penalty scores are doubled at positions 2-13 relative to the 5' end of the small RNA query sequence.
4.  Duplexes are rejected if they:
  * have more than one single-nucleotide bulge or gap.
  * have more than seven total mismatches, G:U base pairs, bulges and gaps.
  * have more than four total mismatches or four total G:U base pairs.

Predicted targets are printed out if they are equal to or lower than the cutoff score specified.

Note: the -i option limits SW to reverse complement matches only, but you can use the -r option with targetfinder.pl to search both strands of a sequence database. This should be done if the database is a genome sequence so that target sites on both strands can be found.

Additional Tools
============
targetfinder_threads.pl
------------
Executes parallel TargetFinder jobs using Perl interpreter threads.

Requirements
------------
1.  Requires Perl 5.10.0 or higher.
2.  Requires the Perl threads and Thread::Queue modules.

Required arguments
------------
  -f <file>      Input small RNA sequences file (FASTA-format).
  
  -d <file>      Target sequence database file (FASTA-format).
  
  -o <file>      Output file. Stores collective results.

Optional Arguments
------------
  -c <float>     Prediction score cutoff value (DEFAULT = 4).
  
  -t <int>       Number of TargetFinder threads/CPUs to use (DEFAULT = 1).
  
  -r             Search reverse strand for targets?. Use this option if the database is genomic DNA.
  
  -h             Shows the help menu.

References
============
1.  Allen E, Xie Z, Gustafson AM, Carrington JC (2005) microRNA-directed phasing during trans-acting siRNA biogenesis in plants. Cell 121: 207Ð221. doi:[10.1016/j.cell.2005.04.004](http://dx.doi.org/10.1016/j.cell.2005.04.004).
2.  Fahlgren N, Howell MD, Kasschau KD, Chapman EJ, Sullivan CM, Cumbie JS, Givan SA, Law TF, Grant SR, Dangl JL, Carrington JC (2007) High-throughput sequencing of *Arabidopsis* microRNAs: Evidence for frequent birth and death of *MIRNA* genes. PLoS ONE 2: e219. doi:[10.1371/journal.pone.0000219](http://dx.doi.org/10.1371/journal.pone.0000219).
3.  Fahlgren N, Carrington JC (2010) miRNA Target Prediction in Plants. Methods in molecular biology (Clifton, NJ) 592: 51Ð57. doi:[10.1007/978-1-60327-005-2_4](http://dx.doi.org/10.1007/978-1-60327-005-2_4).
