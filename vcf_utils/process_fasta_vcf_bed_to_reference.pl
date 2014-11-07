#!/usr/bin/perl
# This script takes a fasta file containing the SNP of interest (flanked by 100bp on either side) and puts in a 
# reference and alternate allele bracket in the middle
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HO7130
use strict;
use Class::Struct;

struct(vcf =>{
	'start' => '$',
	'end' => '$',
	'ref' => '@',
	'alt' => '@',
});

my $usage = "$0 <filtered vcf> <100bp flanking fa> <snp pos (bed) file>\n";
open(VCF, "< $ARGV[0]") || die $usage;
open(FA, "< $ARGV[1]") || die $usage;
open(COM, "< $ARGV[2]") || die $usage;
my @bname = split(/\./, $ARGV[1]);
open(OUT, "> $bname[0]_bracket.fa");

my %lookup; # {chr}->{start} = 1;

# Read in the bed file positions (only the first coordinate and chromosome are important)
while(my $line = <COM>){
	chomp $line;
	$line =~ s/\r//g;
	my @s = split(/\t/, $line);
	$lookup{$s[0]}->{$s[1]} = 1;
}
close COM;

# Now to read in the vcf file and filter out only the regions that we're interested in for analysis
# The bedfile input determines these locations (ie. the %lookup hash)
my %vcfholder; # {chr}->{start} = vcf;
while(my $line = <VCF>){
	chomp $line;
	if($line =~ /^#/){next;}
	my @s = split(/\s+/, $line);
	if($lookup{lc($s[0])}->{$s[1]}){
		my @ref = split(//, $s[3]);
		my @alt = split(//, $s[4]);
		my $start = $s[1] - 100;
		my $end = $s[1] + 200;
		$vcfholder{"chr20"}->{$start} = vcf->new(
			'start' => $start,
			'end' => $end,
			'ref' => \@ref,
			'alt' => \@alt,);
	}
	if(scalar(keys(%{$vcfholder{"chr20"}})) == scalar(keys(%{$lookup{"chr20"}}))){
		last;
	}
}
close VCF;


# Now to take a list of the fasta
my @holder;
my @seq;
while(my $line = <FA>){
	if($holder[0] && $line =~ /^>/){
		my $c = 0;
		my $working = $vcfholder{lc($holder[0])}->{$holder[1]};
		print OUT ">$holder[0]:$holder[1]-$holder[2] $holder[3]\n"; 
		for(my $x = 0; $x < scalar(@seq); $x++){
			my $v = $seq[$x];
			if($c == 100){
				if(!($working->ref->[0] eq $v)){
					print "ERROR for " . join(" ", @holder) . ". REF " . $working->ref->[0] . " does not match base: " . $v . "\n";
				}
				if(scalar(@{$working->ref}) > 1 ){
					$x += scalar(@{$working->ref}) - 1;
				}
				print OUT "[" . join("", @{$working->ref}) . "/" . join("", @{$working->alt}) . "]";
				
			}else{
				print OUT $v;
			}
			$c++;
		}
		print OUT "\n";
		@seq = ();
		chomp $line;
		(@holder) = $line =~ /^>(chr.+):(\d+)-(\d+)(.*)/;
	}elsif($line =~ /^>/){
		chomp $line;
		(@holder) = $line =~ /^>(chr.+):(\d+)-(\d+)(.*)/;
	}else{
		chomp $line;
		$line =~ s/[\[\]]//g;
		push(@seq, split(//, $line));
	}
}

my $c = 0;
my $working = $vcfholder{lc($holder[0])}->{$holder[1]};
print OUT ">$holder[0]:$holder[1]-$holder[2] $holder[3]\n"; 
for(my $x = 0; $x < scalar(@seq); $x++){
	my $v = $seq[$x];
	if($c == 100){
		if(!($working->ref->[0] eq $v)){
			print "ERROR for " . join(" ", @holder) . ". REF " . $working->ref->[0] . " does not match base: " . $v . "\n";
		}
		if(scalar(@{$working->ref}) > 1 ){
			$x += scalar(@{$working->ref}) - 1;
		}
		print OUT "[" . join("", @{$working->ref}) . "/" . join("", @{$working->alt}) . "]";
		
	}else{
		print OUT $v;
	}
	$c++;
}
print OUT "\n";

close FA;
close OUT;
exit;