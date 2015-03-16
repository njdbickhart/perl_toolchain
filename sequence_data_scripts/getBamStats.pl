#!/usr/bin/perl
# This script is designed to run samtools idxstats on each bam input by the user
# The functionality is simple right now, but is designed to be extensible later

use strict;
use Getopt::Std;
use FileHandle;
use bamTools;

my %opts;
my $usage = "perl $0 -g <genome size> -b <bam files (comma separated)> [-o <output tab file [otherwise, STDOUT]>]
REQUIRED:
	-b	Input bam files to check (multiple bams can be input if entered as comma delimited values (no spaces!))
	
OPTIONAL:
	-o	Print data to this output file (output is tab delimited)\n";
	
getopt('bo', \%opts);

if(scalar(@ARGV) == 0 || $ARGV[0] eq '-h'){
	print $usage;
	exit;
}

unless(defined($opts{'b'})){
	print "Missing mandatory arguments!\n";
	print $usage;
	exit;
}

my $fh;
if(defined($opts{'o'})){	
	$fh = FileHandle->new;
	$fh->open("> $opts{o}");
}else{
	$fh = *STDOUT;
}

my $reader = SamFileReader->new();
my $totalCov = 0;
my $bamnum = scalar(split(/,/, $opts{'b'}));
print STDERR "Determining raw x coverage from $bamnum bams...\n";
print $fh "BamName\tRawXCov\tMapXCov\tAvgRawChrcov\tAvgMapChrcov\n";
foreach my $b (split(/,/, $opts{'b'})){
	$reader->prepSam($b);
	my ($rawcov, $mappedcov, $chrcov) = $reader->getXCov();
	$totalCov += $rawcov;
	my ($rawavg, $mappedavg) = averageChrCov($chrcov);
	print $fh "$b\t$rawcov\t$mappedcov\t$rawavg\t$mappedavg\n";
}

if(defined($opts{'o'})){	
	$fh = FileHandle->new;
	$fh->close();
}

exit;

sub averageChrCov{
	my ($hash) = @_;
	my $numchrs = scalar(keys(%{$hash}));
	if($numchrs == 0){
		return 0, 0;
	}
	my $rawsum = 0;
	my $mappedsum = 0;
	foreach my $k (keys(%{$hash})){
		$rawsum += $hash->{$k}->[0];
		$mappedsum += $hash->{$k}->[1];
	}
	$rawsum /= $numchrs;
	$mappedsum /= $numchrs;
	return $rawsum, $mappedsum;
}
