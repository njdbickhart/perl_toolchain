#!/usr/bin/perl
# This script processes a gzipped vcf file to estimate the number of markers per sample, per chromosome

use strict;

chomp(@ARGV);
my $usage = "perl $0 <input gzip vcf> <output tab file>\n";

unless(scalar(@ARGV) == 2){
	print $usage;
	exit;
}

my %data; # {chr}->{sample} = count
my %atab; # unique list of animals
open(my $IN, "gunzip -c $ARGV[0] |");
while(my $line = <$IN>){
	chomp $line;
	if($line =~ /^#CHROM/){
		my @segs = split(/\t/, $line);
		for(my $x = 9; $x < scalar(@segs); $x++){
			$atab{$x} = $segs[$x];
		}
	}elsif($line =~ /^#/){
		next;
	}else{
		my @segs = split(/\t/, $line);
		for(my $x = 9; $x < scalar(@segs); $x++){
			my @gsegs = split(/:/, $segs[$x]);
			unless($gsegs[0] =~ /0\/0/ || $gsegs[0] =~ /\.\/\./){
				$data{$segs[0]}->{$atab{$x}} += 1;
			}
		}
	}
}
close $IN;

open(my $OUT, "< $ARGV[1]");
my @chrs = sort{$a cmp $b} keys(%data);
print {$OUT} "\t" . join("\t", @chrs) . "\n";
foreach my $an (sort{ $a cmp $b} values(%atab)){
	print {$OUT} "$an";
	foreach my $c (@chrs){
		print {$OUT} "\t" . $data{$c}->{$an};
	}
	print {$OUT} "\n";
}
close $OUT;

exit;