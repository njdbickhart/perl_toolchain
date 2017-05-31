#!/usr/bin/perl
# This is a short script designed to output key alignments of sam files without taking up the entire real estate of the screen
# Also, processes alignments to bed file for fast 

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 (-s <input sam> || -b <input bam>) [-d output to bed file]\n";

getopt('sb', \%opts);

unless(defined($opts{'s'}) || defined($opts{'b'})){
	print $usage;
	exit;
}

my $IN;
if(defined($opts{'s'})){
	open($IN, "< $opts{s}") || die "Could not open input sam file: $opts{s}!\n$usage";
}else{
	open($IN, "samtools view $opts{b} |");
}

