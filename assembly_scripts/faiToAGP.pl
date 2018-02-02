#!/usr/bin/perl
# A script that quickly converts a samtools fai file into an AGP for editting
# Useful in conjunction with my fasta editting program
# AGP data:
# contig_15_unplaced      1       27645   1       D       15      49676592        49704237        +
# M       1       16340   1       D       M       1       16340   +
# contig_X_unplaced       1       978900  1       D       contig_X_unplaced       1       978900  +
# Leftover_ScbfJmS_1000   1       41429   1       D       Leftover_ScbfJmS_1000   1       41429   +
# Leftover_ScbfJmS_1002   1       18926   1       D       Leftover_ScbfJmS_1002   1       18926   +

use strict;
my $usage = "perl $0 <fai file> <output agp file>\n";
chomp(@ARGV);
unless(scalar(@ARGV) == 2){
	print $usage;
	exit;
}

open(my $IN, "< $ARGV[0]") || die "Could not open fai file: $ARGV[0]\n$usage";
open(my $OUT, "> $ARGV[1]");
while(my $line = <$IN>){
	chomp $line; 
	my @segs = split(/\t/, $line);
	print {$OUT} "$segs[0]\t1\t$segs[1]\t1\tD\t$segs[0]\t1\t$segs[1]\t+\n";
}
close $IN;
close $OUT;

exit;
