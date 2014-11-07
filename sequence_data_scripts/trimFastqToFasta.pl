#!/usr/bin/perl
# This script takes an input fastq file, calculates the quality score, and trims the beginning and ends of the sequence based on 
# the quality score cutoff input by the user
# This is primarily useful for longer fastqs created by Quiver

use strict;
use Getopt::Std;

my $usage = "perl $0 [options]
	-i	input fastq file
	-s	quality score encoding subtractor (typically \"33\" for Sanger based quality score encodings)
	-t	quality score threshold for trimming ends (after subtraction of the encoding!)
	-o	output fasta file\n";
	
my %opts;
getopt('isto', \%opts);

unless(defined($opts{'i'}) && defined($opts{'s'}) && defined($opts{'t'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

open(IN, "< $opts{i}") || die "Could not open input file!\n$usage";
open(OUT, "> $opts{o}") || die "Could not open output file!\n$usage";
while(my $header = <IN>){
	my $seq = <IN>;
	my $plus = <IN>;
	my $quality = <IN>;
	chomp($seq); chomp($quality); chomp($header);
	
	my ($trimcount, $revised) = sequenceTrimmer($seq, $quality, $opts{'s'}, $opts{'t'});
	$header =~ s/\@/\>/;
	print STDERR "Removed $trimcount bases from: $header\n";
	print OUT "$header\n$revised\n";
}
close IN;
close OUT;

exit;
sub sequenceTrimmer{
	my ($seq, $qual, $offset, $thresh) = @_;
	my @seqs = split(//, $seq);
	my @quals = split(//, $qual);
	my $trimcount = 0;
	
	# Forward direction first
	my $last = 0;
	for(my $x = 0; $x < scalar(@seqs); $x++){
		my $q = ord($quals[$x]) - $offset;
		if($q < 0){
			print STDERR "Error with fastq encoding! $q is less than the offset: $offset!\n";
			exit;
		}
		if($q >= $thresh){
			# this base is greater than the threshold, look ahead three bases to see if we should keep it
			my $keep = 1;
			for(my $y = $x; $y < $x + 3; $y++){
				my $b = ord($quals[$x]) - $offset;
				if($b < $thresh){
					$keep = 0;
					last;
				}
			}
			if($keep){last;}
		}
		$last = $x;
	}
	$trimcount += $last;
	my $num = scalar(@seqs);
	my $newseq = join("", @seqs[$last .. $num]);
	my $newqual = join("", @quals[$last .. $num]);
	
	@seqs = split(//, $newseq);
	@quals = split(//, $newqual);
	
	# Now for the reverse direction
	
	$last = scalar(@seqs);
	for(my $x = scalar(@seqs) - 1; $x > 0; $x--){
		my $q = ord($quals[$x]) - $offset;
		if($q < 0){
			print STDERR "Error with fastq encoding! $q is less than the offset: $offset!\n";
			exit;
		}
		if($q >= $thresh){
			# this base is greater than the threshold, look ahead three bases to see if we should keep it
			my $keep = 1;
			for(my $y = $x; $y > $x - 3; $y--){
				my $b = ord($quals[$x]) - $offset;
				if($b < $thresh){
					$keep = 0;
					last;
				}
			}
			if($keep){last;}
		}
		$last = $x;
	}
	
	$trimcount += scalar(@seqs) - $last;
	return $trimcount, join("", @seqs[0 .. $last]);

}