#!/usr/bin/perl
#########################################################
#  Target Finder
#
#  Copyright 2007 Oregon State University
#
#  Noah Fahlgren
#  Christopher M. Sullivan
#  Kristin D. Kasschau
#  James C. Carrington
#
#  Department of Botany and Plant Pathology
#  Center for Genome Research and Biocomputing
#  Oregon State University
#  Corvallis, OR 97331
#
#  asrp_questions@science.oregonstate.edu
#
# This program is not free software; you can NOT redistribute it and/or
# modify it at all.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#########################################################

=head1 NAME

targetfinder.pl - search for potential miRNA target sites in a sequence database.

=head1 VERSION

See VERSION file.

=head1 REQUIREMENTS

See INSTALL file.

=head1 INSTALL

See INSTALL file.

=head1 INPUT

REQUIRED ARGUMENTS

 -s     Small RNA sequence (either RNA or DNA).
 -d     FASTA-formated sequence file containing potential target sequences.

OPTIONAL ARGUMENTS

 -q     Query sequence name (DEFAULT = query).
 -c     Prediction score cutoff value (DEFAULT = 4).
 -r     Search reverse strand for targets? (BOOLEAN, DEFAULT=FALSE).
 -h     Shows the help menu.

=head1 OUTPUT

targetfinder.pl writes all output to the terminal (STOUT). To save the output to a file use
'>' to redirect output to a file.

ex.

./targetfinder.pl -s UGUGUUCUCAGGUCACCCCUU -d arab_cdna -q miR399a > miR399a_predicted_targets.txt

=head1 CONTACT

Copyright 2007 Oregon State University. Please contact targetfinder@cgrb.oregonstate.edu
for help.

=cut

use strict;
use warnings;
use File::Temp qw(tempfile tempdir);
use Getopt::Std;
use vars qw/ $opt_d $opt_s $opt_q $opt_c $opt_h $opt_r /;
use constant DEBUG => 0;

#########################################################
# Start Variable declarations                           #
#########################################################

# Get environment variables
my $home_dir = $ENV{'HOME'};   # Location of home directory, base directory where temporary directory will be created
#my $fasta34 = $ENV{'FASTA34'}; # Location of the fasta34 binary
my $fasta = "fasta35";

# Creates a temporary directory, unless it exists
my $dir = "$home_dir/tmp";
unless (-e $dir) {
	system "mkdir $home_dir/tmp";  # Note that mkdir must be in your PATH variable
}

if (DEBUG) {
	open (LOG, ">targetfinder.log") or die " Cannot open targetfinder.log: $!\n\n";
}

my ($sRNA, $query_name, $database, $cutoff, $name, @fasta, @fasta_parsed);

getopts('d:s:q:c:hr');
&var_check();

my @tempfileList;

$sRNA =~ tr/atgcu/ATGCU/;

# Short query name for FASTA34 parsing
if ($query_name =~ /(.{1,5})/) {
	$name = $1;
}

print LOG "Small RNA = $sRNA\n" if (DEBUG);
print LOG "Query name = $query_name\n" if (DEBUG);

# miRNA-target alignment scoring matrix
my %bp;
$bp{"AU"} = 0;
$bp{"UA"} = 0;
$bp{"GC"} = 0;
$bp{"CG"} = 0;
$bp{"GU"} = 0.5;
$bp{"UG"} = 0.5;
$bp{"AC"} = 1;
$bp{"CA"} = 1;
$bp{"AG"} = 1;
$bp{"GA"} = 1;
$bp{"UC"} = 1;
$bp{"CU"} = 1;
$bp{"A-"} = 1;
$bp{"U-"} = 1;
$bp{"G-"} = 1;
$bp{"C-"} = 1;
$bp{"-A"} = 1;
$bp{"-U"} = 1;
$bp{"-G"} = 1;
$bp{"-C"} = 1;
$bp{"AA"} = 1;
$bp{"UU"} = 1;
$bp{"CC"} = 1;
$bp{"GG"} = 1;

#########################################################
# End Variable declarations                             #
#########################################################

#########################################################
# Start Main body of Program                            #
#########################################################

# Create temporary input file for FASTA34
my ($input, $infile) = tempfile(DIR=>$dir);
print LOG "Creating temporary file $infile\n" if (DEBUG);
open ($input, ">$infile") or die "Cannot open tempfile for FASTA34 input ($infile): $!\n\n";
push (@tempfileList,$infile);
print $input "\>$query_name\n$sRNA";
close $input;

# Run FASTA
print STDERR " Running $fasta... ";
@fasta = &fasta($infile, $database);
print STDERR "done\n";

# Parse FASTA results
print STDERR " Parsing results... ";
@fasta_parsed = &fasta_parser(@fasta);
print STDERR "done\n";

# Score alignments
print STDERR " Scoring alignments... ";
my @targets = &bp_score(@fasta_parsed);
print STDERR "done\n";

if (@targets) {
	print STDERR " Finding additional hits... ";
	my ($dbase, @additional) = &get_additional($database, @targets);
	while (@additional) {
		push @targets, @additional;
		undef(@additional);
		($dbase, @additional) = &get_additional($dbase, @targets);
	}
	print STDERR "done\n";

	# Get target site coorinates
	print STDERR " Getting target site coordinates... ";
	@targets = get_coords($database, 1, @targets);
	print STDERR "done\n";
}
	
if ($opt_r) {
	print STDERR " Creating reverse database... ";
	my ($rev, $revdb) = tempfile(DIR=>$dir);
	open ($rev, ">$revdb") or die "Cannot open tempfile for reverse database ($revdb): $!\n\n";
	push (@tempfileList,$revdb);
	my ($name, $seq);
	open (DB, $database) or die "Cannot open input database ($database): $!\n\n";
	my $step = 0;
	while (my $line = <DB>) {
		chomp $line;
		if ($step == 0) {
			if ($line =~ /\>(.+)/) {
				$name = $1;
				$step = 1;
			}
		} elsif ($step == 1) {
			if ($line =~ /\>(.+)/ || eof(DB)) {
				if (eof(DB)) {
					$seq .= $line;
				}
				$seq = reverse($seq);
				$seq =~ tr/AGTC/TCAG/;
				print $rev "\>$name\n$seq\n";
				next if (eof(DB));
				$name = $1;
				undef($seq);
			} else {
				$seq .= $line;
			}
		}
	}
	close DB;
	close $rev;
	print STDERR "done\n";
	
	# Run FASTA
	print STDERR " Running FASTA34 on reverse database... ";
	@fasta = &fasta($infile, $revdb);
	print STDERR "done\n";

	# Parse FASTA results
	print STDERR " Parsing reverse results... ";
	@fasta_parsed = &fasta_parser(@fasta);
	print STDERR "done\n";

	# Score alignments
	print STDERR " Scoring reverse alignments... ";
	my @rev_targets = &bp_score(@fasta_parsed);
	print STDERR "done\n";
	if (@rev_targets) {
		print STDERR " Finding additional reverse targets... ";
		my ($dbase, @additional) = &get_additional($revdb, @rev_targets);
		while (@additional) {
			push @rev_targets, @additional;
			undef(@additional);
			($dbase, @additional) = &get_additional($dbase, @rev_targets);
		}
		print STDERR "done\n";
		# Get target site coorinates
		print STDERR " Getting coordinates for reverse targets... ";
		@rev_targets = get_coords($revdb, -1, @rev_targets);
		push @targets, @rev_targets;
		print STDERR "done\n";

	}
}
if (@targets) {
	# Sort output
	@targets = sort by_scores(@targets);
	
	# Print results
	foreach my $line (@targets) {
		print "query=$line->{'query_name'}, target=$line->{'hit_accession'}, score=$line->{'score'}, ";
		print "range=$line->{'target_start'}\-$line->{'target_end'}, strand=$line->{'target_strand'}\n\n";
		print "target  5' $line->{'target_seq'} 3'\n";
		print "           $line->{'homology_string'}\n";
		print "query   3' $line->{'miR_seq'} 5'\n\n";
	}
} else {
	print "No results for $query_name\n";
}

#remove all of the temporary files
foreach my $fname (@tempfileList){
	unlink "$fname";
}
exit;

#########################################################
# End Main body of Program                              #
#########################################################

#########################################################
# Start Subroutines                                     #
#########################################################

# Run FASTA
sub fasta {
	print LOG "Running $fasta...\n" if (DEBUG);
	my $input = shift;
	my $db = shift;
	my @output;
#	open FASTA, "$fasta -n -H -Q -f -16 -r +15/-10 -g -10 -w 100 -W 25 -b 100000 -i -U $input $db 1 |";
	open FASTA, "$fasta -n -H -Q -f -16 -r +15/-10 -g -10 -w 100 -W 25 -E 100000 -i -U $input $db 1 |";
	while (<FASTA>) {
		print LOG $_ if (DEBUG);
		push (@output, $_);
	}
	close FASTA;
	return @output;
}

# Parse FASTA output
sub fasta_parser {
	print LOG "Parse results...\n" if (DEBUG);
	my @input = @_;
	my $hit_accession;
	my $id;
	my $spacer;
	my $miRNA;
	my $count = 0;
	my $length = 0;
	my $target;
	my $step = 0;
	my @output;
	foreach my $line (@input) {
		chomp $line;
		if ($step == 0) {
			if ($line =~ /^\>{2}(.+)\(\d+ nt\)/) {
				$hit_accession = $1;
				if ($hit_accession =~ /(^.{6})/) {
					$id = $1;
				}
				while ($hit_accession =~ /\s$/) {
					$hit_accession =~ s/\s$//g;
				}
				$step = 1;
				print LOG "     HIT ACCESSION = $hit_accession\n" if (DEBUG);
			}
		} elsif ($step == 1) {
			if ($line =~ /($name\-\s+)\b(\D+)\b/) {
				$spacer = $1;
				$miRNA = $2;
				$count = ($spacer =~ tr/[A-Z][a-z][0-9][\-][\ ][_]//);
				$miRNA =~ s/\s//g;
				$length = length $miRNA;
				$miRNA =~ tr/AUGC/UACG/;
				print LOG "     miRNA = $miRNA\n" if (DEBUG);
			} elsif ($line =~ /$id\s+/) {
				if ($line =~ /.{$count}(\D{$length})/) {
					$target = $1;
					$target =~ s/\s//g;
					print LOG "     TARGET = $target\n" if (DEBUG);
				}
				push @output, "$hit_accession\t$miRNA\t$target";
				$step = 0;
			}
		}

	}
	return @output;
}

# Score FASTA alignments
sub bp_score {
	print LOG "Score alignments...\n" if (DEBUG);
	my @input = @_;
	my $counter = 0;
	my @output;
	foreach my $line (@input) {
		chomp $line;
		my $mismatch = 0;
		my $gu = 0;
		my $total_mispair = 0;
		my ($hit_accession, $miR_seq, $target_seq) = split /\t/, $line;
		my $miR_length = length $miR_seq;
		my $target_length = length $target_seq;
		next if ($target_seq !~ /^[AGCU-]+$/);
		next if ($miR_length != $target_length);
		next if ($miR_seq =~ /\-.{0,}\-/);
		next if ($target_seq =~ /\-.{0,}\-/);
		next if ($miR_seq =~ /\-/ && $target_seq =~ /\-/);
		my @miRNA = split //, $miR_seq;
		my @target = split //, $target_seq;
		my $cycle = 0;
		my $score = 0;
		my $old_score = 0;
		my $homology_string;
		for (1..$miR_length) {
			$cycle++;
			my $miR_base = pop @miRNA;
			my $target_base = pop @target;
			if ($cycle == 1) {
				my $position = $bp{"$miR_base$target_base"};
				if ($position == 1) {
					$mismatch++;
					$homology_string .= ' ';
				} elsif ($position == 0.5) {
					$gu++;
					$homology_string .= '.';
				} else {
					$homology_string .= ':';
				}
				$score = $position;
			} elsif ($cycle > 13) {
				my $position = $bp{"$miR_base$target_base"};
				if ($position == 1) {
					$mismatch++;
					$homology_string .= ' ';
				} elsif ($position == 0.5) {
					$gu++;
					$homology_string .= '.';
				} else {
					$homology_string .= ':';
				}
				my $new_score = $position;
				$old_score = $score;
				$score = ($old_score + $new_score);
			} else {
				my $position = ($bp{"$miR_base$target_base"}*2);
				if ($position == 2) {
					$mismatch++;
					$homology_string .= ' ';
				} elsif ($position == 1) {
					$gu++;
					$homology_string .= '.';
				} else {
					$homology_string .= ':';
				}
				my $new_score = $position;
				$old_score = $score;
				$score = ($old_score + $new_score);
			}


		}
		$total_mispair = ($mismatch+$gu);
		next if ($total_mispair > 7 || $mismatch > 4 || $gu > 4);
		if ($score <= $cutoff) {
			my %hash;
			$counter++;
			$hash{'hit_accession'} = $hit_accession;
			$hash{'miR_seq'} = $miR_seq;
			$hash{'query_name'} = $query_name;
			$hash{'target_seq'} = $target_seq;
			$hash{'score'} = $score;
			$hash{'homology_string'} = reverse $homology_string;
			$hash{'target_start'} = 0;
			$hash{'target_end'} = 0;
			$hash{'target_strand'} = 0;
			push @output, \%hash;
		}
	}
	return @output;
}

# Look for additional target sites in identified targets
sub get_additional {
	my ($db, @found) = @_;
	my ($fh, $filename) = tempfile(DIR=>$dir);
	open ($fh, ">$filename") or die "Cannot open tempfile for additional hits database ($filename): $!\n\n";
	push (@tempfileList,$filename);
	my $step = 0;
	my $id;
	my $seq;
	open (DB, $db) or die "Cannot open previous database ($db): $!\n\n";
	while (my $line = <DB>) {
		chomp $line;
		if ($step == 0) {
			if ($line =~ /\>(.+)/) {
				$id = $1;
				$step = 1;
			}
		} elsif ($step == 1) {
			if ($line =~ /\>(.+)/ || eof(DB)) {
				if (eof(DB)) {
					if ($line !~ /^\s*$/) {
						$seq .= $line;
					}
				}
				my $match = 0;
				foreach my $hit (@found) {
					my $info = $hit->{'hit_accession'};
					$info =~ s/\\/\\\\/g;
					$info =~ s/\|/\\\|/g;
					$info =~ s/\(/\\\(/g;
					$info =~ s/\)/\\\)/g;
					$info =~ s/\[/\\\[/g;
					$info =~ s/\]/\\\]/g;
					$info =~ s/\{/\\\{/g;
					$info =~ s/\}/\\\}/g;
					$info =~ s/\^/\\\^/g;
					$info =~ s/\$/\\\$/g;
					$info =~ s/\*/\\\*/g;
					$info =~ s/\+/\\\+/g;
					$info =~ s/\?/\\\?/g;
					$info =~ s/\./\\\./g;
					if ($id =~ /$info/) {
						$match++;
						my $hit_seq = $hit->{'target_seq'};
						$hit_seq =~ s/U/T/g;
						$hit_seq =~ s/-//g;
						my $length = length($hit_seq);
						my $mask;
						for (1..$length) {
							$mask .= 'N';
						}
						#$seq =~ s/$hit_seq/$mask/;
						$seq =~ s/$hit_seq/$mask/i;
					}
				}
				if ($match > 0) {
					print $fh "\>$id\n$seq\n";
				}
				$id = $1;
				undef($seq);
			} else {
				$seq .= $line;
			}
		}
	}
	close $fh;
	@fasta = &fasta($infile, $filename);
	@fasta_parsed = &fasta_parser(@fasta);
	my @results = &bp_score(@fasta_parsed);
	#if ($db ne $database) {
	#	unlink $db;
	#}
	if (!@results) {
		unlink $filename;
	} else {
		$db = $filename;
	}
	return ($db, @results);
}

sub get_coords {
	my ($db, $strand, @found) = @_;
	my (%db, $name, $seq);
	open (DB, $db) or die "Cannot open sequence database ($db): $!\n\n";
	my $step = 0;
	while (my $line = <DB>) {
		chomp $line;
		if ($step == 0) {
			#if ($line =~ /\>(.+)/) {
			if ($line =~ /\>(.{1,94})/) {
				$name = $1;
				while ($name =~ /\s$/) {
					$name =~ s/\s$//g;
				}
				$step = 1;
			}
		} elsif ($step == 1) {
			#if ($line =~ /\>(.+)/ || eof(DB)) {
			if ($line =~ /\>(.{1,94})/ || eof(DB)) {
				if (eof(DB)) {
					$seq .= $line;
				}
				$db{$name} = $seq;
				next if (eof(DB));
				$name = $1;
				while ($name =~ /\s$/) {
					$name =~ s/\s$//g;
				}
				undef($seq);
			} else {
				$seq .= $line;
			}
		}
	}
	close DB;
	foreach my $hit (@found) {
		my $seq = $db{$hit->{'hit_accession'}};
		my $hit_seq = $hit->{'target_seq'};
		$hit_seq =~ s/\-//g;
		$hit_seq =~ s/U/T/g;
		my $length = length($hit_seq);
		my $offset;
		if ($strand == 1) {
			if ($seq =~ /^(.*)$hit_seq/i) {
				$offset = length($1);
				$hit->{'target_start'} = (length($1) + 1);
				$hit->{'target_end'} = ($hit->{'target_start'} + $length - 1);
			}
		} elsif ($strand == -1) {
			if ($seq =~ /(.*)$hit_seq(.*)$/) {
				$offset = length($1);
				$hit->{'target_start'} = (length($2) + 1);
				$hit->{'target_end'} = ($hit->{'target_start'} + $length - 1);
			}
		}
		my $mask;
		for (1..$length) {
			$mask .= 'N';
		}
		substr($seq,$offset,$length,$mask);
		$db{$hit->{'hit_accession'}} = $seq;
		$hit->{'target_strand'} = $strand;
	}
	return @found;
}

# Sort output by score, acsending
sub by_scores {
	$a->{'score'} <=> $b->{'score'}
}

# Check options
sub var_check {
	if ($opt_s) {
		$sRNA = $opt_s;
	} else {
		&var_error();
	}
	if ($opt_d) {
		$database = $opt_d;
	} else {
		&var_error();
	}
	if ($opt_q) {
		$query_name = $opt_q;
	} else {
		$query_name = 'query';
	}
	if ($opt_c) {
		$cutoff = $opt_c;
	} else {
		$cutoff = 4;
	}
	if ($opt_h) {
		&var_error();
	}
	if (!$home_dir) {
		print " ERROR: environmental variable 'HOME' not set\n";
		&var_error();
	}
	if (!defined($fasta)) {
		print " ERROR: environmental variable 'FASTA' not set\n";
		&var_error();
	}
}

# Print help
sub var_error {
	print "\n\n";
	print " targetfinder.pl will predict potential miRNA target sites using methods described in Allen et al. 2005 and Fahlgren et al. 2007.\n";
	print " You did not provide enough information\n";
	print " Usage: targetfinder.pl -s <sequence> -d <target database> [OPTIONS]\n";
	print " REQUIRED:\n";
	print " -s     Small RNA sequence\n";
	print " -d     Target database path (FASTA-format)\n";
	print " OPTIONAL:\n";
	print " -q     Query sequence name (DEFAULT = query)\n";
	print " -c     Score cutoff value (DEFAULT = 4)\n";
	print " -r     Search reverse strand for targets? (BOOLEAN, DEFAULT=FALSE)\n";
	print " -h     Show this menu\n\n";
	print " Type perldoc targetfinder.pl for more help.\n";
	print "\n\n";
	exit 0;
}

#########################################################
# End Subroutines                                       #
#########################################################
