#!/usr/bin/perl
# This script takes a vcf file and a samtools indexed fasta file and generates a file for use by GeneSeek for probe construction
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HO7130
use strict;
use Class::Struct;
use Getopt::Std;

select((select(STDOUT), $| = 1)[0]);

struct(vcf =>{
	'start' => '$',
	'end' => '$',
	'ref' => '@',
	'alt' => '@',
	'qual' => '$',
	'seq' => '$',
});

my %opts;
getopt('vfoin', \%opts);
my $usage = "perl $0 -v <vcf file> -i <samtools indexed fasta> -o <output fasta file name> -f <optional filter phrase> -n <optional number start>\n";

unless(defined($opts{'v'}) && defined($opts{'i'}) && defined($opts{'o'})){
	print $usage;
	exit;
}
my $filter = 0;
my $fname;
if(defined($opts{'f'})){
	$filter = 1;
	$fname = $opts{'f'}; # just to avoid confusing Perl when I use this in a regex
}

# I'm going to have to run through the VCF twice, otherwise I'll hit penalties for memory and time
# First time, select the variants that I want to save
my %vcfholder; # {chr}->{start} = vcf;
my $found = 0;
my $lnum = 0;
open(VCF, "< $opts{v}") || die "Could not open VCF File: $opts{v}!\n";
while(my $line = <VCF>){
	chomp $line;
	if($line =~ /^#/){next;}
	my @s = split(/\s+/, $line);
	if($s[7] =~ /^INDEL.+/){next;} # avoid INDELs like the plague;
	$lnum++;
	if($lnum % 100000 == 0){
		print "Processed up to line:\t$lnum\tfound: $found\r";
	}
	if($filter){
		if(!($line =~ /$fname/)){
			next;
		}else{
			$found += 1;
		}
	}else{
		$found += 1;
	}
	
	my @ref = split(//, $s[3]);
	my @alt = split(//, $s[4]);
	my $start = $s[1] - 101;
	my $end = $s[1] + 100;
	$vcfholder{$s[0]}->{$s[1]} = vcf->new(
		'start' => $start,
		'end' => $end,
		'ref' => \@ref,
		'alt' => \@alt,
		'qual' => $s[5]);
	$vcfholder{$s[0]}->{$s[1]}->seq(samtools_faidx_wrapper($s[0], $vcfholder{$s[0]}->{$s[1]}, $opts{'i'}));
	
}
seek(VCF, 0, 0);
print "\n";

# Now we begin to traverse the VCF file to find the common variant sites and annotate them
my @store; #I'm going to do a "lookback analysis" to try to catch the target variants in this array.
my $lastchr = "hey";
$lnum = 0;
while(my $line = <VCF>){
	my @segs = split(/\t/, $line);
	if($line =~ /^#/){next;}
	if($segs[7] =~ /^INDEL.+/){next;} # avoid INDELs like the plague;
	$lnum++;
	
	if($lnum % 100000 == 0){
		print "Processed up to line:\t$lnum\r";
	}
	if($segs[0] ne $lastchr){
		if(scalar(@store) > 1){
			@store = ();
		}
		$lastchr = $segs[0];
	}
	if(exists($vcfholder{$segs[0]}->{$segs[1]}) && scalar(@store) < 25){
		# In this case, we don't have 25 variants so we have to prematurely process this one
		my $num = scalar(@store) - 1;
		for(my $x = 0; $x < 25; $x++){
			my $l = <VCF>;
			my @s = split(/\t/, $line);
			my $temp = convert_vcf_line_to_object(\@s);
			push(@store, $temp);
		}
		$vcfholder{$segs[0]}->{$segs[1]} = modify_fasta_seq(\@store, $vcfholder{$segs[0]}->{$segs[1]}, $segs[1]);
	}elsif(scalar(@store) > 25){
		my $vstart = $store[25]->start();
		if(exists($vcfholder{$segs[0]}->{$vstart})){
			$vcfholder{$segs[0]}->{$vstart} = modify_fasta_seq(\@store, $vcfholder{$segs[0]}->{$vstart}, $vstart);
		}
	}
	push(@store, convert_vcf_line_to_object(\@segs));
	if(scalar(@store) >= 50){
		# We only want to keep 50 vcf lines in memory at a time
		shift(@store);
	}
}
print "\n";
close VCF;

# Now to print out the SNPs
print "Printing SNPs\n";
my $number = 0;
if(defined($opts{'n'})){
	$number = $opts{'n'};
}
open(OUT, "> $opts{o}");
foreach my $chr (sort{$a cmp $b} keys(%vcfholder)){
	foreach my $s (sort {$a <=> $b} keys(%{$vcfholder{$chr}})){
		if($vcfholder{$chr}->{$s}->qual() < 99){next;} # Low quality site filter
		my $seq = $vcfholder{$chr}->{$s}->seq();
		my $ns = ($seq =~ tr/NX/NX/);
		if($ns > 25){next;} # Too many repetitive/anonymous bases
		
		$number++;
		my $numstr = sprintf("%06d", $number);
		my $fname = "ARS-USDA-AGIL-$chr\-$s\-$numstr";
		print OUT ">$fname\n";
		print OUT "$seq\n";
	}
}
close OUT;

exit;
sub vcf_object_tostr{
	my ($vcf) = @_;
	my $str = $vcf->start() . "\t";
	$str .= $vcf->end() . "\t";
	$str .= $vcf->seq() . "\t";
	$str .= join(";", @{$vcf->ref()}) . "\t";
	$str .= join(";", @{$vcf->alt()}) . "\t";
	return $str;
}
# This subroutine will modify the nucleotide sequence surrounding the SNP
# It will also tag the SNP site by adding brackets around it an displaying the ref and alt alleles
sub modify_fasta_seq{
	my ($store, $variant, $avoid) = @_;
	my $start = $variant->start();
	my $end = $variant->end();
	my $seq = $variant->seq();
	
	# Let's convert all of the variant locations within the bounds of the probe sequences
	foreach my $s (@{$store}){
		my $vstart = $s->start();
		if($vstart != $avoid && $vstart >= $start && $vstart <= $end){
			my $pos = $vstart - $start;
			if($pos > length($seq)){
				print "\nAttempted substr of $seq at position: $pos. $vstart - $start\n";
				print vcf_object_tostr($variant);
				print "\n";
			}
			my $cbase = nucleotide_conversion_tool($s->ref->[0], @{$s->alt()});
			substr($seq, $pos, 1) = $cbase;
		}
	}
	
	# Now to highlight the SNP site
	my $refb = $variant->ref->[0];
	my $altb = $variant->alt->[0];
	my $absbase = $avoid - $start;
	substr($seq, $absbase, 1) = "[$refb\\$altb]";
	$variant->seq($seq);
	return $variant;
}

# Simple sub to try to quickly convert variants into objects
sub convert_vcf_line_to_object{
	my ($aref) = @_;
	my @ref = split(//, $aref->[3]);
	my @alt = split(//, $aref->[4]);
	return vcf->new(
		'start' => $aref->[1],
		'ref' => \@ref,
		'alt' => \@alt,
		'qual' => $aref->[5]);
}
# This wrapper gets the samtools fasta sequence for the vcf object
sub samtools_faidx_wrapper{
	my ($chr, $vcf, $fasta) = @_;
	my $s = $vcf->start();
	my $e = $vcf->end();
	open(SAM, "samtools faidx $fasta $chr:$s\-$e |");
	my $h = <SAM>;
	my $seq;
	while(my $line = <SAM>){
		chomp $line;
		$seq .= $line;
	}
	close SAM;
	return $seq;
}

# This conversion tool should help me mark known variant sites within the 100 bp of the probe
sub nucleotide_conversion_tool{
	my %iupac = (
		"AG" => "R",
		"CT" => "Y",
		"CG" => "S",
		"AT" => "W",
		"GT" => "K",
		"AC" => "M",
		"CGT" => "B",
		"AGT" => "D",
		"ACT" => "H",
		"ACG" => "V");
	my ($ref, $alt, $salt) = @_;
	my @base = ($ref, $alt);
	if(defined($salt)){
		push(@base, $salt);
	}
	@base = sort{$a cmp $b} @base;
	my $code = join("", @base);
	if(exists($iupac{$code})){
		return $iupac{$code};
	}else{
		return "N";
	}
}