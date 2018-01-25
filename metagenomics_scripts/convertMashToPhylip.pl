#!/usr/bin/perl
# This script is designed to take a folder of ".dist" files from MASH and consolidate them into a Phylip-style matrix

use strict;
use Getopt::Std;

my $usage = "perl $0 -f <input folder of dist files> -o <output phylip formatted distance matrix>
NOTE: files must have the same PATH in the header of the dist matrices.\n";

my %opts;
getopt('fo', \%opts);

unless(defined($opts{'f'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my @dists = `ls $opts{f}/*.dist`;
if(scalar(@dists) < 1){
	print "Did not find any distance files!\n$usage";
	exit;
}

my @dataOrder;
#my %idLookup;
chomp(@dists);
# Get header and generate unique IDs for each sample
open(my $IN, "< $dists[0]") || die "Could not open initial distance file!\n";
my $head = <$IN>;
chomp $head;
my @hsegs = split(/\s+/, $head); 
foreach my $h (@hsegs){
	if($h eq "#query"){next;}
	my $id = GetUniqueID($h);
	push(@dataOrder, $id);
}
close $IN;

my %data;
# now to process each file and concatenate the data
foreach my $d (@dists){
	open(my $IN, "< $d") || die "Could not open distance file: $d\n";
	my $h = <$IN>;
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\s+/, $line);
		my $id = GetUniqueID($segs[0]);
		$segs[0] = $id;
		my @holder;
		for(my $x = 1; $x < scalar(@segs); $x++){
			my $f = sprintf("%.4f", $segs[$x]);
			push(@holder, $f . "");
		}
		$data{$id} = \@holder;
	}
	close $IN;
}

# Finally, print out the Phylip output
open(my $OUT, "> $opts{o}");
my $num = scalar(@dataOrder);
print {$OUT} "$num\n";
foreach my $j (@dataOrder){
	my $rows = $data{$j};
	my $twhite = join(" ", @{$rows});
	my $id = sprintf("%-10s", $j);
	print {$OUT} "$id $twhite\n";
}
close $OUT;

exit;

sub GetUniqueID{
	my ($name) = @_;
	my @nsegs = split(/\//, $name);
	my ($fname) = $nsegs[-1] =~ /(^.+)\..+$/; # remove file extension
	if(length($fname) > 10){
		my $num = length($fname);
		my $temp = substr($fname, length($fname) - 10);
		$fname = $temp;
	}	
	return $fname;
}
