#!/usr/bin/perl
# This script is designed to run samtools idxstats on each bam input by the user
# The functionality is simple right now, but is designed to be extensible later

use strict;
use Getopt::Std;
use FileHandle;
use bamTools;

my %opts;
my $usage = "perl $0 -b <bam files (comma separated)> || -l <ls of bam files> || -n <newline delimited list of bams>[-o <output tab file [otherwise, STDOUT]>]
REQUIRED:
	-b	Input bam files to check (multiple bams can be input if entered as comma delimited values (no spaces!))
	
OPTIONAL:
	-o	Print data to this output file (output is tab delimited)\n";
	
if(scalar(@ARGV) == 0 || $ARGV[0] eq '-h'){
	print $usage;
	exit;
}

getopt('boln', \%opts);



unless(defined($opts{'b'}) || defined($opts{'l'}) || defined($opts{'n'})){
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

my @bams;
if(defined($opts{'b'})){
	@bams = split(/,/, $opts{'b'});
}elsif(defined($opts{'l'})){
	@bams = `ls $opts{l}`;
	if(scalar(@bams) < 1){
		print STDERR "Error finding bam files from the following LS command: $opts{l}!\n";
		exit;
	}
}elsif(defined($opts{'n'})){
	open(my $IN, "< $opts{n}") || die "Could not open list of delimited bams!\n$usage";
	while(my $line = <$IN>){
		chomp $line;
		push(@bams, $line);
	}
	close $IN;
}
print STDERR "Determining raw x coverage from $bamnum bams...\n";
print $fh "BamName\tTotalReads\tMappedReads\tUnmappedReads\tRawXCov\tMapXCov\tAvgRawChrcov\tAvgMapChrcov\n";
foreach my $b (@bams){
	$reader->prepSam($b);
	my ($mappedReads, $unmappedReads) = $reader->getReads();
	my ($rawcov, $mappedcov, $chrcov) = $reader->getXCov();
	$totalCov += $rawcov;
	my ($rawavg, $mappedavg) = averageChrCov($chrcov);
	my $totalReads = $mappedReads + $unmappedReads;
	print $fh "$b\t$totalReads\t$mappedReads\t$unmappedReads\t$rawcov\t$mappedcov\t$rawavg\t$mappedavg\n";
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
