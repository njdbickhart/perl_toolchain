#!/usr/bin/perl
# This is your first and final warning -- don't use this for data to be published!
# If you do, you will be cursed to a thousand years of programming errors. You will always forget semicolons and your variable names will be
# indecipherable! Trust me, it will be terrible!

use strict;
my $usage = "perl $0 <first otu table> <second otu table> <first shared file> <second shared file> <output basename>\n";
chomp @ARGV;
unless(scalar(@ARGV) == 5){
	print $usage;
	exit;
}

# Read in OTU table one, tally counts
open(my $IN, "< $ARGV[0]");
my $head = <$IN>;
my %otu1;
my %otuc1;
my %otuc2;
my %otuF;
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$otu1{$segs[2]} = $segs[1];
	$otuc1{$segs[2]} = $segs[0];
}
close $IN;

# Read in OTU table 2, tally final counts
open(my $IN, "< $ARGV[1]");
my $head = <$IN>;
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$otuF{$segs[2]} += $segs[1];
	$otuc2{$segs[2]} = $segs[0];
}
close $IN;

# Rerank OTUs and rescore
my @sortedOTU = sort{$otuF{$b} <=> $otuF{$a}} keys(%otuF);
my %conv1;
my %conv2;
my @oturanks;
my %otutax;
for(my $x = 0; $x < scalar(@sortedOTU); $x++){
	my $otunum = $x + 1;
	my $oname = "Otu$otunum";
	$otuF{$oname} = $otuF{$sortedOTU[$x]};
	$otutax{$oname} = $sortedOTU[$x];
	$conv1{$otuc1{$sortedOTU[$x]}} = $oname;
	$conv2{$otuc2{$sortedOTU[$x]}} = $oname;
	push(@oturanks, $oname);
}

# Now read in shared file and print out final shared file
open(my $OUT, "> $ARGV[4].shared");
print {$OUT} "label\tGroup\tnumOtus\t" . join('\t', @oturanks) . "\n";
open(my $IN, "< $ARGV[2]");
my $head = <$IN>;
my @hsegs = split(/\t/, $head);
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	print {$OUT} "$segs[0]\t$segs[1]\t$segs[2]\t";
	my @retOtus = getOtuCounts(\@oturanks, \%conv1, \@hsegs, \@segs);
	print {$OUT} join("\t", @retOtus) . "\n";
}

close $IN;

open(my $IN, "< $ARGV[3]");
my $head = <$IN>;
my @hsegs = split(/\t/, $head);
while(my $line = <$IN>){
	chomp $line;
		my @segs = split(/\t/, $line);
		print {$OUT} "$segs[0]\t$segs[1]\t$segs[2]\t";
		my @retOtus = getOtuCounts(\@oturanks, \%conv2, \@hsegs, \@segs);
		print {$OUT} join("\t", @retOtus) . "\n";
}
close $IN;
close $OUT;

# Now update the final taxonomy table
open(my $OUT, "> $ARGV[4].taxonomy");
print {$OUT} "OTU\tSize\tTaxonomy\n";
foreach my $o (@oturanks){
	print {$OUT} "$o\t" . $otuF{$o} . "\t" . $otutax{$o} . "\n";
}
close $OUT;

exit;

sub getOtuCounts{
	my ($ranks, $conv, $hsegs, $segs) = @_;
	my %oIndex;
	for(my $x = 3; $x < scalar(@{$hsegs}); $x++){
		$oIndex{$conv->{$hsegs->[$x]}} = $x;
	}
	my @otus;
	foreach my $r (@{$ranks}){
		if(!exists($oIndex{$r})){
			push(@otus, 0);
		}else{
			push(@otus, $segs->[$oIndex{$r}]);
		}
	}
	return @otus;
}