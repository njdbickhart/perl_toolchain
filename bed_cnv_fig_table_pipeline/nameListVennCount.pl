#!/usr/bin/perl
# This script is designed to take two files, identify the unique names, and then peform venn counts
# I'm going to start off easy at first
# Adding the ability to print out groups to separate files for later retrieval

use strict;

chomp(@ARGV);

my $usage = "Must enter at least two files\nperl $0 [-o : flag to print out groupings to text] file1 file2 file3 ... fileN\n";

my $outputflag = 0;
if($ARGV[0] eq "-o"){
	$outputflag = 1;
	shift(@ARGV);
}

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
my %printout;
foreach my $k (keys(%store)){
	my $aref = $store{$k};
	my %set;
	foreach my $a (@{$aref}){
		$set{$a} = 1;
	}
	
	my $group = join(";", sort{ $a <=> $b} keys(%set));
	if($outputflag){
		push(@{$printout{$group}}, $k);
	}
	$counter{$group} += 1;
}

print "Set\tCount\n";
foreach my $k (sort { $a cmp $b} keys(%counter)){
	print "$k\t$counter{$k}\n";
}

if($outputflag){
	print "\n";
	foreach my $g (keys(%printout)){
		my $n = $g;
		$n =~ s/;/_/g;
		open(OUT, "> group_$n.txt");
		print OUT join("\n", @{$printout{$g}});
		print OUT "\n";
		close OUT;
		print "Group: $g in output file: group_$n.txt\n";
	}
}

exit;