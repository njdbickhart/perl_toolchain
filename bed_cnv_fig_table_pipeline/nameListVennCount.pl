#!/usr/bin/perl
# This script is designed to take two files, identify the unique names, and then peform venn counts
# I'm going to start off easy at first

use strict;

chomp(@ARGV);

my $usage = "Must enter at least two files\nperl $0 file1 file2 file3 ... fileN\n";

if(scalar(@ARGV) <= 1){
	print $usage;
	exit;
}

my %store;
for(my $x = 0; $x < scalar(@ARGV); $x++){
	my $fnum = $x + 1;
	print "File Number $fnum: $ARGV[$x]\n";
	open(IN, "< $ARGV[$x]") || die "Could not open file $fnum : $ARGV[$x]!\n";
	while(my $line = <IN>){
		chomp $line;
		push(@{$store{$line}}, $fnum);
	}
	close IN;
}

my %counter;
foreach my $k (keys(%store)){
	my $aref = $store{$k};
	my %set;
	foreach my $a (@{$aref}){
		$set{$a} = 1;
	}
	
	my $group = join(";", sort{ $a <=> $b} keys(%set));
	$counter{$group} += 1;
}

print "Set\tCount\n";
foreach my $k (sort { $a cmp $b} keys(%counter)){
	print "$k\t$counter{$k}\n";
}

exit;