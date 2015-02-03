#!/usr/bin/perl
# This script will sample the bam file originally used to generate RAPTR-SV files and use it to filter RAPTR-SV calls
# Filtration is based on an assumed X% coverage of the original bam (lowest coverage for assumed heterozygote calls 
# based on simple Poisson coverage estimates) and a X% support ratio (support decimal divided by total number of reads)

use strict;
use Getopt::Std;
use bamTools;
use RaptrSVUtils;

my %opts;
my $usage = "perl $0 -b <sam/bam file (list?)> -i <RAPTR-SV file> -o <output bed file> -f <support ratio [optional]> -c <coverage filter [optional]> 
	-b	REQUIRED: sam/bam files used to generate this RAPTR-SV output file (can be a comma delimited list: ie. bam1,bam2,bam3)
	-i	REQUIRED: input; the RAPTR-SV output file
	-o	REQUIRED: the name of the output bed file
	
	-f	OPTIONAL: Read support ratio [default: 0.50]
	-c	OPTIONAL: X coverage filter [default: calls must have at least 0.30 of the expected coverage in read pairs]\n";
	
getopt('bifco', \%opts);

if(!defined($opts{'b'}) || !defined($opts{'i'}) || !defined($opts{'o'})){
	print "Error! Missing required arguments!\n"; 
	print $usage;
	exit;
}

my $supFilt = (defined($opts{'f'}))? $opts{'f'} : 0.50;
my $covFilt = (defined($opts{'c'}))? $opts{'c'} : 0.30;

# First, sample the coverage of the bams and determine the coverage threshold
my $reader = SamFileReader->new();
my $totalCov = 0;
my $bamnum = scalar(split(/,/, $opts{'b'}));
print STDERR "Determining raw x coverage from $bamnum bams...\n";
foreach my $b (split(/,/, $opts{'b'})){
	my ($rawcov, $mappedcov) = $reader->getXCov($b);
	$totalCov += $rawcov;
	print STDERR "BAM $b had $rawcov Raw X coverage and $mappedcov Mapped X coverage\n";
}

my $thresh = $totalCov * $covFilt;
my $filteredEvents = 0;
my $totalLines = 0;
print STDERR "Read coverage threshold: $thresh\n";
print STDERR "Now filtering the RAPTR-SV file...\n";
open(IN, "< $opts{i}") || die "Could not open RAPTR-SV file: $opts{i}!\n$usage";
open(OUT, "> $opts{o}");
while(my $line = <IN>){
	$totalLines++;
	my $sv = SVEntry->new('line' => $line);
	my $rawcount = $sv->getRawReadCount();
	my $ratio = $sv->getSupportRatio();
	
	if($rawcount >= $thresh && $ratio >= $supFilt){
		print OUT $sv->printBed . "\n";
	}else{
		$filteredEvents++;
	}
}
close IN;
close OUT;
print STDERR "Filtered $filteredEvents out of $totalLines total events\n";

exit;
