#!/usr/bin/perl
# This script is designed to pull sam file strings from a bam file
# Very useful if you want to recover read pairs from a bam when you know the mappings of one
# of the mates

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -i <input bam> -r <input, newline separated read names> -o <output sam>\n";

getopt('iro', \%opts);

unless(defined($opts{'i'}) && defined($opts{'r'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# Place reads into a hash for fast access
my %reads;
open(IN, "< $opts{r}") || die "Could not open read name file!\n";
while(my $line = <IN>){
	chomp $line;
	$reads{$line} = 1;
}
close IN;

# Now, just read the bam and pull the sam strings we're looking for
open(IN, "samtools view $opts{i} |") || die "Could not read bam file!\n";
open(OUT, "> $opts{o}");
while(my $line = <IN>){
	# Not chomping as that might take more time for a long file
	my @segs = split(/\t/, $line);
	if(exists($reads{$segs[0]})){
		print OUT $line;
	}
}
close IN;
close OUT;
print STDERR "Finished scanning bam\n";
exit;