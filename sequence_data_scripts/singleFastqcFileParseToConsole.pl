#!/usr/bin/perl
# This script will process an individual fastq file to test my fastqcParser module
# Output is to STDOUT, tab-delimited

use strict;
use simpleLogger;
use fastqcParser;

my $usage = "perl $0 <fastqc exe> <fastq file> ... <fastq file, N>  : pushes all output to standard output\n";

chomp(@ARGV); 

unless(scalar(@ARGV) > 1){
	print $usage;
	exit;
}

my $log = simpleLogger->new('logFileBaseStr' => 'simple_parser');

my @parsers;
my $fqcexe = shift(@ARGV);
foreach my $f (@ARGV){
	#file => fastq file, sample => sample_name, library => library_name, readNum => first_or_second_read, log => simpleLogger
	my $fqcparse = fastqcParser->new('file' => $f, 'sample' => 'sample', 'library' => 'lib', 'readNum' => '1', 'log' => $log);
	$fqcparse->runFastqc($fqcexe);
	$fqcparse->parseStats();
	push(@parsers, $fqcparse);
}

my @headers = $parsers[0]->getHeaderArray();
print join("\t", @headers) . "\n";
foreach my $p (@parsers){
	my @output = $p->getOutArray();
	print join("\t", @output) . "\n";
}

exit;