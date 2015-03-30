#!/usr/bin/perl
# This script processes a list of preprocessed bams associated with animal names in order to call
# individual VCF segments and create intersection files

use strict;
use Getopt::Long;
use seqHaplotypeTools;
use bamTools;

# required input files
my ($animaltable, $haplotypetable, $bamlist, $reffasta); 

# method flags
my $snpfilter = 0;
my $snpcompare = 0;
my $verbose = 0;

my $usage = "perl $0 --animal <full animal table> --haplotype <haplotype coords> --bams <bam table> --fasta <reference fasta>

	--animal	Tab delimited table containing selected animal data
	--haplotype	Tab delimited table containing haplotype locations in genome
	--bams		Tab delimited file with file loc, followed by inferred animal name (from animal data table)
	--fasta		Reference genome fasta file
	
	OPTIONAL:
	--filtersnps	Flag: use perl-based tool to filter (place in \'.filter\' file) heterozygous SNPs
	--comparesnps	Flag: use VCFtools to compare SNPs from different animals (prior to filtering)
	--verbose	Flag: debugging or extra information\n";

GetOptions("animal=s" => \$animaltable,
	"haplotype=s" => \$haplotypetable,
	"bams=s" => \$bamlist,
	"fasta=s" => \$reffasta,
	"comparesnps" => \$snpcompare,
	"filtersnps" => \$snpfilter,
	"verbose" => \$verbose)
	|| die "Missing mandatory arguments!\n$usage";
	

foreach my $prog (@{"samtools", "bwa"}){
	StaticUtils::checkReqs($prog);
}

if(!defined($animaltable) && !defined($haplotypetable) && !defined($bamlist)){
	print "Missing required input data tables!\n";
	print $usage;
	exit;
}

# Load necessary input tables, animals.info and haplotype.info
my $hapSegs = HapSegMap->new();
$hapSegs->loadFromFile($haplotypetable);

my $anData = Individual->new();
$anData->loadFromFile($animaltable);

print STDERR "[HAPCALLER] Finished loading data tables!\n";

# Generate hash of list of animal bams
my %bamFiles;
open(IN, "< $bamlist") || die "[HAPCALLER] error opening bam file table!\n$usage";
while(my $line = <IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$bamFiles{$segs[1]} = $segs[0];
}
close IN;

# Now, calculate which pairings we need to make
my @presentAns = keys(%bamFiles);
$anData->calcIntersections();
my($multiple, $single) = $anData->returnIntersections(\@presentAns);

print STDERR "[HAPCALLER] Of " . scalar(@presentAns) . " total animals, found " . scalar(@{$single}) . " singleton haps and " . scalar(@{$multiple}) . "\n";

# OK, now the engine will start processing the files to generate VCFs
# Single first
my $samtools = SamtoolsExecutable->new();
foreach my $sing (@{$single}){
	# get hap coords
	my $hap = $sing->hap->hap;
	my $allele = $sing->hap->allele;
	my $hapseg = $sing->getS($hap);
	my $ucsc = $hapseg->printUCSC();
	my $output = $sing->animals->[0] . ".$hap.$allele.vcf";
	
	$samtools->GenerateSamtoolsVCF($sing->animals, $output, $reffasta, $ucsc);
}

my @multFiles;
foreach my $mult (@{$multiple}){
	# get hap coords
	my $hap = $sing->hap->hap;
	my $allele = $sing->hap->allele;
	my $hapseg = $sing->getS($hap);
	my $ucsc = $hapseg->printUCSC();
	my $tempAnStr = join("-", @{$mult-animals});
	my $output = "$tempAnStr.$hap.$allele.vcf";
	
	push(@multFiles, $output);
	
	$samtools->GenerateSamtoolsVCF($mult->animals, $output, $reffasta, $ucsc);
}

if($snpfilter){
	print STDERR "[HAPCALLER] SNP filter not implemented yet!\n";
}

if($snpcompare){
	print STDERR "[HAPCALLER] SNP comparison not implemented yet!\n";
}

exit;