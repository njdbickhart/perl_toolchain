#!/usr/bin/perl
# Base algorithm adapted from: http://genomics-array.blogspot.com/2011/02/calculating-n50-of-contig-assembly-file.html

use strict;

my $usage = "perl $0 <input fasta file>\n";
chomp(@ARGV);
unless(scalar(@ARGV) == 1){
	print $usage;
	exit;
}

open(my $IN, "< $ARGV[0]") || die "Could not open input fasta file: $ARGV[0]!\n";

## Read Fasta File and compute length ###
my $length;
my $totalLength; 
my @arr;
while(my $line = <$IN>){
	chomp $line; 
	if($line =~ /^>/){
		push (@arr, $length);
		$totalLength += $length; 
		$length = 0;
		next;
	}
	$length += length($line);
}

close($IN);

my @sort = sort {$b <=> $a} @arr;
my $n50; 
my $L50;
foreach my $val(@sort){
	$n50+=$val;
	$L50++;
	if($n50 >= $totalLength/2){
		print "N50 length:\t$n50\n";
		print "N50 value:\t$val\n";
		print "L50 value:\t$L50\n";
		last; 
	}
}

exit;
