#!/usr/bin/perl
# This script is designed to get the sam stats and list them in a counter table for rapid diagnostics

use strict;

our %flags = (
	0 => 'TotalReads',
	1 => 'MultSegs',
	2 => 'BothSegsAlignProperly',
	4 => 'SegUnmapped',
	8 => 'SecondSegUnmapped',
	16 => 'ThisSeqReverse',
	32 => 'NextSeqReverse',
	64 => 'ThisIsFirstRead',
	128 => 'ThisIsSecondRead',
	256 => 'SecondaryAlignment',
	512 => 'FailedQCFilter',
	1024 => 'Duplicate',
	2048 => 'SuplementaryAlign' );

my @nums = (0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048);	
my $usage = "perl $0 <Bam1> <Bam2> ... <BamN>\n";
chomp(@ARGV);

unless(scalar(@ARGV) > 0){
	print $usage;
	exit;
}


foreach my $b (@ARGV){
	print "File: $b\n\n";
	
	my %table;
	open(IN, "samtools view $b |") || die "Could not open bam file: $b!\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		foreach my $k (@nums){
			if(($segs[1] & $k) == $k){
				$table{$k} += 1;
			}
		}
	}
	close IN;
	printMrkDownTable(\%table);
}

exit;

sub printMrkDownTable{
	my ($hash) = @_;
	my $nlen = 25;
	
	my $vlen = 0;
	foreach my $v (values(%{$hash})){
		if(length($v) > $vlen){
			$vlen = length($v);
		}
	}
	$vlen++;
	
	printf "\|%-*s\|%*s\|\n", $nlen, "Flag", $vlen, "Count";
	my $sepstr = '-' x ($nlen - 1);
	my $sepcon = '-' x ($vlen - 1);
	printf("\|\:%s\|%s\:\|\n", $sepstr, $sepcon);
	
	foreach my $k (sort {$a <=> $b} keys(%{$hash})){
		my $ftranslate = $flags{$k};
		my $count = $hash->{$k};
		printf "\|%-*s\|%*d\|\n", $nlen, $ftranslate, $vlen, $count;
	}
	print "\n";
}
