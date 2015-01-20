#!/usr/bin/perl
# This script converts AIPL genotypes into ped map format files
# There is an option ("a") to use a list of alt/ref allele codes instead of A A, and B B notation
# Output from this option requires a bed file with the ref/alt codes in the 4th column, unseparated (eg. "AC")

use strict;
use snpGenotypeAPI;
use Getopt::Std;

my %opts;
my $usage = "perl $0 [options]
	-c <chromosome.data file> 
	-g <AIP genotype file> 
	-o <output base filename> 
	-p <population string> 
	-x <gender assignment (can be a number (1=male, 2=female) or a tab delimited file (anid\tgendercode) <OPTIONAL>>
	-a <use actual ref/alt alleles in this file <OPTIONAL>>
	-n <prefilter chromosome.data with chip name <OPTIONAL>>\n";

getopt('cgopanx', \%opts);
my $population = (defined($opts{'p'}))? $opts{'p'} : "NA";
my $codes = (defined($opts{'a'}))? $opts{'a'} : "NO";

unless(defined($opts{'c'}) && defined($opts{'g'}) && defined($opts{'o'})){
	print $usage;
	exit;
}
my $gender = (defined($opts{'x'}))? $opts{'x'} : "other";

my $aiplmap = mapBed->new();
if(defined($opts{'n'})){
	$aiplmap->convertFromCHRDATA($opts{'c'}, $opts{'n'});
}else{
	$aiplmap->convertFromCHRDATA($opts{'c'});
}
print STDERR "Loaded chromosome.data file...\n";

my @genotypes;

if($codes ne "NO"){
	print STDERR "Creating Genotype code!\n";
	my $coder = CreateGenotypeCode->new();
	$coder->createGenoCodeMap($codes, $aiplmap, "codes.geno");
}

open(IN, "< $opts{g}") || die "Could not open AIP genotype file!\n";
while(my $line = <IN>){
	chomp $line;
	my $temp = AIPLGenotype->new();
	$temp->convertGenoStr($line, $population);
	print STDERR "Working on animal: " . $temp->ankey() . "\n";
	$temp->maplocs($aiplmap);
	push(@genotypes, $temp);
}
close IN;

my $printer = PedMapPrinter->new();
if($codes ne "NO"){
	$printer->printOutPedMap(\@genotypes, $opts{'o'}, $gender, "codes.geno");
}else{
	$printer->printOutPedMap(\@genotypes, $opts{'o'}, $gender);
}
print STDERR "Finished creating output. Output is in files: $opts{o}.ped and $opts{o}.map\n";
exit;