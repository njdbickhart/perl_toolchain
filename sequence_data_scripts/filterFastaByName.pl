#!/usr/bin/perl
# This is a simple utility designed to filter away fasta entries by name

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -f <input fasta file> -l <new-line delimited list of names to filter> -o <output fasta file>\n";

getopt('flo', \%opts);

if(!defined($opts{'f'}) || !defined($opts{'l'})){
	print $usage;
	exit;
}

my %remove;
open(my $IN, "< $opts{l}") || die "Could not open list of names to filter!\n$usage";
while(my $line = <$IN>){
	chomp $line;
	$line =~ s/>//g;
	$remove{$line} = 1;
}
close $IN;

my $filter = 0;
my $numfiltered = 0;
open(my $OUT, "> $opts{o}");
open(my $IN, "< $opts{f}") || die "Could not open fasta file to filter!\n$usage";
while(my $line = <$IN>){
	chomp $line;
	
	if($line =~ /^>/){
		$line =~ s/>//g;
		if(exists($remove{$line})){
			$filter = 1;
			$numfiltered++;
		}else{
			print $OUT ">$line\n";
			$filter = 0;
		}
	}elsif($filter){
		next;
	}else{
		print $OUT "$line\n";
	}
}
close $IN;
close $OUT;

print STDERR "Filtered $numfiltered Fasta entries!\n";

exit;