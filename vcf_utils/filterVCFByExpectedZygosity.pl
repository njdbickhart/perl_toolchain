#!/usr/bin/perl
# This script is designed to take input samples and only select variants that are
# consistent with expected genotype status. Values that allow for "fuzzy" genotypes
# are also used.

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 [options]
	-v	input vcf file to filter
	-p	input individual file. Names match chrom line. Tab delimited \"het\", \"hom\", \"wild\"
	-o	output vcf file name
	-u	UCSC region to focus on in output file (optional)
	-f	fuzzy heterozygous search -- will count homozygote alts as well (optional)
";

getopt('vpou', \%opts);

if(!defined($opts{'v'}) || !defined($opts{'p'}) || !defined($opts{'o'})){
	print $usage;
	exit;
}

my $regionFilter = (defined($opts{'u'}))? 1 : 0;
my ($chr, $start, $end);
if($regionFilter){
	($chr, $start, $end) = $opts{'u'} =~ m/(.+):(\d+)-(\d+)/;
}

# Load zygosity records
my %samples;
open(my $IN, "< $opts{p}") || die "Could not open individual file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$samples{$segs[0]} = $segs[1];
}
close $IN;

# Now to open the vcf and figure out which columns to ID
my %cols;
my @search;
open($IN, "< $opts{v}") || die "Could not open vcf file!\n";
open(my $OUT, "> $opts{o}");
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[0] eq "#CHROM"){
		for(my $x = 9; $x < scalar(@segs); $x++){
			if(exists($samples{$segs[$x]})){
				$cols{$x} = $samples{$segs[$x]};
				push(@search, $x);
			}
		}
		print $OUT "$line\n";
		last;
	}
	print $OUT "$line\n";
}

# Now to plough through the variant calls
my $threshold = scalar(@search);
my $found = 0;
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($regionFilter){
		if($segs[0] ne $chr && ($segs[1] < $start || $segs[1] > $end)){
			next; # skip this if it's outside of our region of interest
		}
	}
	
	my $obs = 0;
	foreach my $ind (@search){
		my $exp = $cols{$ind};
		my @vcfsegs = split(/:/, $segs[$ind]);
		
		if($exp eq "het"){
			if($vcfsegs[0] =~ /1[\/|]0/ || $vcfsegs[0] =~ /0[\/|]1/
				|| ($opts{'f'} && $vcfsegs[0] =~ /1[\/|]1/)){
				$obs++;
			}
		}elsif($exp eq "hom"){
			if($vcfsegs[0] =~ /1[\/|]1/){
				$obs++;
			}
		}else{
			# Wildtype expected
			if($vcfsegs[0] =~ /0[\/|]0/){
				$obs++;
			}
		}
	}
	if($obs >= $threshold){
		print $OUT "$line\n";
		$found++;
	}
}
close $IN;
close $OUT;

print STDERR "Identified $found lines that fit expected zygosity profiles\n";

exit;