#!/usr/bin/perl
# This is a script designed to process a tab file with sequence data fastqs and produce slurm scripts

use strict;
use Getopt::Std;
use slurmTools;

my %opts;
my $usage = "perl $0 -b <base outfolder name> -o <base output file name> -t <input tab sequence files>\n";
getopt('bot', \%opts);

unless(defined($opts{'b'}) && defined($opts{'o'}) && defined($opts{'t'})){
	print $usage;
	exit;
}


