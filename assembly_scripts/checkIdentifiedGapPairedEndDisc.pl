#!/usr/bin/perl
# This script processes the tab delimited output from identifyFilledGaps.pl, a bam file containing paired end read alignments to the target assembly
# and generates a score for gap filling based on soft clipping and discordant read mappings
# OUTPUT:
# Class gap_name gap_length close_coordinates close_len read_num softclipped_bases oae_reads
#
# class can be: "FullClose, CrypticMis, GAP, Unknown"

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -g <tab delimited gap closure file> -t <target reference bam file> -o <output tab file summary>\n";

getopt('gto', \%opts);

unless(defined($opts{'g'}) && defined($opts{'t'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

unless( -s $opts{'t'} && -s "$opts{t}.bai"){
	print "Error! Bam file either doesn't exist or isn't indexed by samtools!\n";
	exit;
}

my $tested = 0; my $skipped = 0;

open(my $IN, "< $opts{g}") || die "Could not open gap closure file!\n";
open(my $OUT, "> $opts{o}");
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my $gapname = "$segs[1]:$segs[2]-$segs[3]";
	my $gaplen = $segs[3] - $segs[2];
	my ($alignc1, $aligns1, $aligne1) = convert_ucsc($segs[6]);
	my ($alignc2, $aligns2, $aligne2) = convert_ucsc($segs[7]);
	if($segs[0] eq "Closed"){
		my ($tstart, $tend) = get_coords($aligns1, $aligns2, $aligne1, $aligne2);
		my $closecoords = "$segs[5]:$tstart-$tend";
		my $closelen = $tend - $tstart;
		my $type = "Unknown";
		if($segs[8] == -1){
			# We have a very close overlap closure
			# We need to test this to see if its correct or a misassembly			
			# padding start and end coords with 5 bases to check around alignment
			my ($reads, $sclip, $oae) = process_bam_file($segs[5], $tstart - 5, $tend + 5);
			if($reads > 1 && ($sclip > 1 || $oae > 1)){
				$type = "CrypticMis";
			}elsif($reads > 1){
				$type = "FullClose";
			}
			
			print {$OUT} "$type\t$gapname\t$gaplen\t$closecoords\t$closelen\t$reads\t$sclip\t$oae\n";
		}elsif(abs($segs[8] - $segs[4]) > 100000){
			# The difference in closure size is too great here. There's an alignment issue or misassembly
			print {$OUT} "$type\t$gapname\t$gaplen\t$closecoords\t$closelen\t-1\t-1\t-1\n";
		}elsif($segs[8] != $segs[10]){
			# It's pretty clear, this is an unclosed gap!
			$type = "GAP";
			print {$OUT} "$type\t$gapname\t$gaplen\t$closecoords\t$closelen\t-1\t-1\t-1\n";
		}else{
			# This should handle all reasonable length closure events
			my ($reads, $sclip, $oae) = process_bam_file($segs[5], $tstart, $tend);
			if($reads > 1 && ($sclip > 1 || $oae > 1)){
				$type = "CrypticMis";
			}elsif($reads > 1){
				$type = "FullClose";
			}
			print {$OUT} "$type\t$gapname\t$gaplen\t$closecoords\t$closelen\t$reads\t$sclip\t$oae\n";
		}
		$tested++;
	}else{
		$skipped++;
	}
	if($tested % 100 == 0){
		print ".";
	}
}
print "\n";
close $IN;
close $OUT;

print "Tested: $tested\n";
print "Skipped: $skipped\n";
		
exit;
sub process_bam_file{
	my ($bam, $chr, $start, $end) = @_;
	my $reads = 0; my $sclip = 0; my $oae = 0;
	open(my $BAM, "samtools view $bam $chr:$start-$end |");
	while(my $sam = <$BAM>){
		chomp $sam;
		my @segs = split(/\t/, $sam);
		# check if this is a oae
		if($segs[6] eq "*"){
			$oae++;
		}
		
		# process cigar looking for a count of softclipped bases
		while($segs[5] =~ m/(\d+)(\D{1})/g){
			my $cnum = $1;
			my $cstr = $2;
			if($cstr eq "H" || $cstr eq "S"){
				$sclip += $cnum;
			}
		}
		$reads++;
	}
	close $BAM;
	return $reads, $sclip, $oae;
}

sub convert_ucsc{
	my ($ucsc) = @_;
	my @segs = split(/[:-]/, $ucsc);
	return $segs[0], $segs[1], $segs[2];
}

sub get_coords{
	my (@positions) = @_;
	@positions = sort {$a <=> $b} @positions;
	return $positions[1], $positions[2];
}