#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use threads;
use Thread::Queue;
use FindBin qw($Bin);

################################################################################
# Begin variable declarations
################################################################################

my (%opt, $file, $database, $threads, $cutoff, $outfile, $format);

getopts('f:t:d:c:o:p:rh', \%opt);
var_check();

################################################################################
# End variable declarations
################################################################################

################################################################################
# Begin main program
################################################################################

# Load sequences
our %rna : shared;
my $previous;
open(IN, $file) or die "Cannot open file $file: $!\n\n";
while (my $line = <IN>) {
	chomp $line;
	next if ($line =~ /^\s*$/);
	if (substr($line,0,1) eq '>') {
		$previous = substr($line,1);
		$rna{$previous} = '';
	} else {
		$rna{$previous} .= $line;
	}
}
close IN;

# Initialize and run threads
my $workq = Thread::Queue->new();
my $retq = Thread::Queue->new();
while (my ($name, $seq) = each(%rna)) {
	my $job = "$Bin/targetfinder.pl -s $seq -d $database -q $name -c $cutoff -p $format";
	if ($opt{'r'}) {
		$job .= ' -r';
	}
	$workq->enqueue($job);
}
$workq->enqueue("EXIT") for (1 .. $threads);
threads->create("targetfinder") for (1 .. $threads);

# Wait for threads to finish
open(OUT, ">$outfile") or die "Cannot open $outfile: $!.\n\n";
print OUT "{\n" if ($format eq 'json');
while (threads->list(threads::running)) {
	while ($retq->pending) {
		my $result = $retq->dequeue;
		if ($format eq 'json') {
			my @lines = split /\n/, $result;
			shift @lines;
			pop @lines;
			print OUT join("\n", @lines)."\n";
		} else {
			print OUT $result;
		}
	}
	sleep 10;
}
while ($retq->pending) {
	my $result = $retq->dequeue;
	if ($format eq 'json') {
		my @lines = split /\n/, $result;
		shift @lines;
		pop @lines;
		print OUT join("\n", @lines);
	} else {
		print OUT $result;	
	}
}
print OUT "\n}\n" if ($format eq 'json');
close OUT;

################################################################################
# End main program
################################################################################

################################################################################
# Begin subroutines
################################################################################

########################################
# Function: targetfinder
#     Execute TargetFinder threads
########################################
sub targetfinder {
	while (my $job = $workq->dequeue()) {
		last if $job eq 'EXIT';
		my $result = `$job`;
		$retq->enqueue($result);
	}
	threads->detach;
}

########################################
# Function: var_check
#     Direct input arguments to the
#     correct variables and check for
#     missing arguments
########################################
sub var_check {
	if ($opt{'h'}) {
		var_error();
	}
	if ($opt{'f'}) {
		$file = $opt{'f'};
	} else {
		var_error();
	}
	if ($opt{'d'}) {
		$database = $opt{'d'};
	} else {
		var_error();
	}
	if ($opt{'c'}) {
		$cutoff = $opt{'c'};
	} else {
		$cutoff = 4;
	}
	if ($opt{'t'}) {
		$threads = $opt{'t'};
	} else {
		$threads = 1;
	}
	if ($opt{'o'}) {
		$outfile = $opt{'o'};
	} else {
		var_error();
	}
	if ($opt{'p'}) {
		$format = $opt{'p'};
		unless ($format eq 'classic' || $format eq 'gff' || $format eq 'json' || $format eq 'table') {
			print STDERR "Output format $format is not valid.\n\n";
			exit 1;
		}
	} else {
		$format = 'classic';
	}
}

########################################
# Function: var_error
#     Error message/instructions to print
########################################
sub var_error {
	print STDERR "\n\n";
	print STDERR "TargetFinder (threads): Plant small RNA target prediction tool, parallelized with threads.\n\n";
	print STDERR "Usage:   targetfinder_threads.pl -f <FASTA file> -d <target database> -o <output file> [options]\n\n";
	print STDERR "Options: -f <file>    Input small RNA sequences file (FASTA-format)\n";
	print STDERR "         -d <file>    Target sequence database file (FASTA-format)\n";
	print STDERR "         -o <file>    Output file. Stores collective results\n";
	print STDERR "         -c <float>   Prediction score cutoff value (DEFAULT = 4)\n";
	print STDERR "         -t <int>     Number of TargetFinder threads/CPUs to use (DEFAULT = 1)\n";
	print STDERR "         -p <str>     Output format for small RNA-target pairs (DEFAULT = 'classic')\n";
	print STDERR "                      Available options: 'classic' (Original TargetFinder base-pairing format)\n";
	print STDERR "                                         'gff'     (Generic Feature Format)\n";
	print STDERR "                                         'json'    (JavaScript Object Notation)\n";
	print STDERR "                                         'table'   (Tab-deliminated Format)\n";
	print STDERR "         -r           Search reverse strand for targets?. Use this option if the database is genomic DNA.\n";
	print STDERR "         -h           Print this menu\n";
	print STDERR "\n\n";
	exit 1;
}

################################################################################
# End subroutines
################################################################################
