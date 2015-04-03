#!/usr/bin/perl
# This script is designed to work on up to 1000 bed files at once, with the same genome coordinates
# It opens all of these bed files and reads each line to see if there is a uniform lack of coverage
# (or lower coverage) within the regions

use strict;
use Getopt::Std;
use FileHandle;

my $usage = "perl $0 -i <comma separated list of bams> -t <integer coverage threshold> -o <output>\n";
my %opts;

getopt('ito', \%opts);

unless(defined($opts{'i'}) && defined($opts{'t'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $thresh = $opts{'t'};
my @files = split(/,/, $opts{'i'});

my @FH;
foreach my $f (@files){
	my $fh = FileHandle->new();
	$fh->open("< $f");
	push(@FH, $fh);
}

open(OUT, "> $opts{o}");
my $num = scalar(@FH);
my $startfh = $FH[0];
while(my $line = <$startfh>){
	my @segs = split(/\t/, $line);
	my $chr = $segs[0];
	my $start = $segs[1];
	my $end = $segs[2];
	my $vUnderThresh = 0;
	
	if($segs[-1] <= $thresh){
		$vUnderThresh++;
	}
	for(my $x = 1; $x < $num; $x++){
		my $fh = $FH[$x];
		$line = <$fh>;
		@segs = split(/\t/, $line);
		if($segs[-1] <= $thresh){
			$vUnderThresh++;
		}
	}
	if($vUnderThresh == $num){
		print OUT "$chr\t$start\t$end\n";
	}
}
close OUT;

exit;
