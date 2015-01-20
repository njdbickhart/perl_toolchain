#!/usr/bin/perl
# This is a one-shot script designed to use my snp genotype API to convert AIPl formats over to ped and map files

use strict;
use snpGenotypeAPI;

my $population = "BTHO";

my $infile = "genotypes.dat";
my $inmap = "/home/dbickhart/share/umd3_data/snplocs/chromosome.data";
my $fiftyk = "/home/dbickhart/share/umd3_data/snplocs/illumina_50k_umd3.bed";
my $hd = "/home/dbickhart/share/umd3_data/snplocs/illumina_HD_umd3.bed";
my $genocode = "overlap.geno";
my $topalleles = "/home/dbickhart/share/umd3_data/snplocs/aiplallelecode.bed";
my $outbase = "dprholstein";

my $aiplmap = mapBed->new();
$aiplmap->convertFromCHRDATA($inmap);

my $fiftymap = mapBed->new();
$fiftymap->convertFromBed($fiftyk);

my $hdmap = mapBed->new();
$hdmap->convertFromBed($hd);

my @genotypes;
my $createdcode = 0;
open(IN, "< $infile") || die "Could not open $infile!\n";
while(my $line = <IN>){
	chomp $line;
	my $temp = AIPLGenotype->new();
	$temp->convertGenoStr($line, $population);
	print STDERR "Working on animal: " . $temp->anname() . "\n";
	my ($newMap, $newgen) = $aiplmap->intersectLocs($fiftymap, $temp->geno());
	my ($finalMap, $finalGen) = $newMap->intersectLocs($hdmap, $newgen);
	
	if(!$createdcode){
		print STDERR "Creating Genotype code!\n";
		my $coder = CreateGenotypeCode->new();
		$coder->createGenoCodeMap($topalleles, $finalMap, $genocode);
		$createdcode = 1;
	}
	
	$temp->maplocs($finalMap);
	$temp->geno($finalGen);
	push(@genotypes, $temp);
}
close IN;

my $printer = PedMapPrinter->new();
$printer->printOutPedMap(\@genotypes, $outbase, $genocode);
exit;