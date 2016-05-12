#!/usr/bin/perl
# This is a small front-end script to quickly process a VCF file into tabular output

use strict;
use Getopt::Std;
use vcfTools;

my $usage = "perl $0 -f <gzipped, indexed VCF file> -o <output tab file> -u <[OPT] UCSC coordinates to search> -a <[OPT] new-line delimited list file of animals to include>\n";

my %opts;
getopt('foua', \%opts);

unless(defined($opts{'f'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $worker = VCFFilterAndSplit->new('vcffile' => $opts{'f'}, 'output' => $opts{'o'});
# Add optional components
if(defined($opts{'u'})){
	$worker->ucsc($opts{'u'});
}
if(defined($opts{'a'})){
	$worker->indivs($opts{'a'});
}

# Run the main routine
$worker->generateSummaryTabFile();

print STDERR "Program output in: $opts{o}\n";
exit;