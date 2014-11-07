#!/usr/bin/perl
# This script is designed to take only the lines that have SNPS in the input animals (derived from the header section of the VCF)
# Works only on the output from the GATK at the moment
# Does not remove other calls at this point (though it should be trivial to do that in the future)

use strict;
my $usage = "Usage: $0 <input file> <animal1> <animal2> ... <output vcf>
	[note: can use \"all\" to denote the selection of all snps instead of an animal range]\n";

chomp(@ARGV);
if(scalar(@ARGV) < 3){
	print $usage;
	exit;
}

my $input = shift(@ARGV);
my $output = pop(@ARGV);
my @indivs = @ARGV;

open(IN, "gunzip -c $input |") || die "Could not open input: $input!\n $usage";
open(OUT, "> $output") || die "Could not open output: $output!\n $usage";

my @index;
while(my $line = <IN>){
	if($line =~ /^##/){
		print OUT $line;
	}elsif($line =~ /^#CHROM/){
		chomp $line;
		my @segs = split(/\s+/, $line);
		for(my $x = 9; $x < scalar(@segs); $x++){
			foreach my $v (@indivs){
				if($segs[$x] eq $v || $v eq "all"){
					push(@index, $x);
					last;
				}
			}
		}
		if(scalar(@index) == 0){
			die "Could not find any input animals in the VCF header columns!\n";
		}
		print OUT "$line\n";
	}else{
		if(scalar(@index) == 0){
			die "Could not find any input animals in the VCF header columns or could not find the VCF header!\n";
		}
		chomp $line;
		my @segs = split(/\s+/, $line);
		my $print = 0;
		foreach my $i (@index){
			my @s = split(/:/, $segs[$i]);
			if($s[0] eq "1/1" || $s[0] eq "0/1"){
				$print = 1;
				last;
			}
		}
		if($print){
			print OUT "$line\n";
		}
	}
}

close IN;
close OUT;
exit;