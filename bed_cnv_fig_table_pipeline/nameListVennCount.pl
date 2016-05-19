#!/usr/bin/perl
# This script is designed to take two files, identify the unique names, and then peform venn counts
# I'm going to start off easy at first
# Adding the ability to print out groups to separate files for later retrieval
# 5/19/2016:	v2	Added ability to print a specific grouping to STDOUT

use strict;

chomp(@ARGV);

my $usage = "$0 v2\nMust enter at least two files\nperl $0 [-o OR -l <groups, comma delimited>] file1 file2 file3 ... fileN
	-o	[Optional] Must be the first argument. Prints output to separate files.
	-l	[Optional] Must be the first argument. Prints output from a specific group (underscore notation based on file order) to STDERR\n";

my $outputflag = 0;
my $stdoutflag = 0;
my %stdoutgroup;
my %involvedFiles;
if($ARGV[0] eq "-o"){
	$outputflag = 1;
	shift(@ARGV);
}elsif($ARGV[0] eq "-l"){
	$stdoutflag = 1;
	shift(@ARGV);
	my $groups = shift(@ARGV);
	my @subgroups = split(/,/, $groups);
	foreach my $s (@subgroups){
		# Order the numbers in proper format
		my @bs = split(/_/, $s);
		# Add file numbers that are involved in the estimation
		foreach my $j (@bs){
			$involvedFiles{($j + 1)} = 1;
		}
		if(scalar(@bs) == 1){
			$stdoutgroup{$bs[0]} = 1;
		}else{
			my $sg = join(";", @bs);
			$stdoutgroup{$sg} = 1;
		}
	}
}

if(scalar(@ARGV) <= 1){
	print $usage;
	exit;
}

my %store;
for(my $x = 0; $x < scalar(@ARGV); $x++){
	my $fnum = $x + 1;
	unless($stdoutflag && !exists($involvedFiles{$fnum})){
		print "File Number $fnum: $ARGV[$x]\n";
	}
	open(IN, "< $ARGV[$x]") || die "Could not open file $fnum : $ARGV[$x]!\n";
	while(my $line = <IN>){
		chomp $line;
		$line =~ s/\r//g;
		# Store each unique element in a file in an array with file number associated with it
		# NOTE: redundant entries in a file may be counted twice!
		push(@{$store{$line}}, $fnum);
	}
	close IN;
}

my %counter;
my %printout;
my $stdoutcount = 0;
#my %uniqueList;
foreach my $k (keys(%store)){
	my $aref = $store{$k};
	my %set;
	foreach my $a (@{$aref}){
		$set{$a} = 1;
	}
	
	my $group = join(";", sort{ $a <=> $b} keys(%set));
	if($outputflag){
		# We want to print these out in the end, so store them in a list to print later
		push(@{$printout{$group}}, $k);
	}elsif($stdoutflag && exists($stdoutgroup{$group})){
		# We just want the output from this list, so print the data now
		print STDERR "$k\n";
		$stdoutcount++;
	}
	$counter{$group} += 1;
}

if($stdoutflag){
	# If we didn't print anything, notify the user
	if($stdoutcount == 0){
		print "EMPTY\n";
	}
	exit;
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