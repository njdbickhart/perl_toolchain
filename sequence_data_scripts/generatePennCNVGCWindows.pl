#!/usr/bin/perl
# This script is designed to take mapped SNP probes and generate GC window values around them

use strict;
use Getopt::Std;
use fastaTools;

my %opts;
my $usage = "perl $0 -s <snp list (blast format)> -r <reference genome fasta> -w <window radius around SNP (bp)> -o <output file>\n";
getopt('srwo', \%opts);

unless(defined($opts{'s'}) && defined($opts{'r'}) && defined($opts{'w'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# initiate fastaTools
my $worker = fastaTools->new('fasta' => $opts{'r'});

# load SNP list and calculate GC content for each one
open(my $OUT, "> $opts{o}");
open(my $IN, "< $opts{s}") || die "Could not open SNP list!\n$usage";
my $h = <$IN>;
print {$OUT} "Name\tChr\tPosition\tGC\n";
while(my $line = <$IN>){
	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line);
	my $start = ($segs[4] - $opts{'w'} < 0)? 1 : $segs[4] - $opts{'w'};
	my $end = $segs[4] + $opts{'w'};
	my $gc = $worker->calcGCContent($segs[2], $start, $end);
	
	$gc *= 100;
	
	print {$OUT} "$segs[1]\t$segs[2]\t$segs[4]\t$gc\n";
}
close $IN:
close $OUT;

exit;