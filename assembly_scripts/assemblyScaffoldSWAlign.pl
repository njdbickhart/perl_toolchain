#!/usr/bin/perl
# This script takes contig names from several assemblies and attempts to align them using a 
# modified Smith-Waterman alignment algorithm. No Affine penalty.

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -a <comparison file one> -b <comparison file two> -o <output alignment>
	The comparison files must be an ordered list of names, separated by newlines
	The B file is reversed to provide a reverse alignment\n";
	
getopt('abo', \%opts);

unless(defined($opts{'a'}) && defined($opts{'b'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# Read in files
my @adata;
my @bdata;

open(my $IN, "< $opts{a}");
while(my $line = <$IN>){
	chomp $line;
	push(@adata, $line);
}
close $IN;

open($IN, "< $opts{b}");
while(my $line = <$IN>){
	chomp $line;
	push(@bdata, $line);
}
close $IN;

my @forwardmatrix = CreateSWMatrix(\@adata, \@bdata);

my ($xpos, $ypos, $highscore) = FindMaxValue(\@forwardmatrix);

my (@foraalign, @forbalign);
Alignment(\@forwardmatrix, $xpos, $ypos, \@adata, \@bdata, \@foraalign, \@forbalign);

open(my $OUT, "> $opts{o}");
PrintAlignments(\@foraalign, \@forbalign, $OUT, 0);

my @revb = reverse(\@bdata);

my @reversematrix = CreateSWMatrix(\@adata, \@revb);
($xpos, $ypos, $highscore) = FindMaxValue(\@reversematrix);
my (@revaalign, @revbalign);
Alignment(\@reversematrix, $xpos, $ypos, \@adata, \@revb, \@revaalign, \@revbalign);
PrintAlignments(\@revaalign, \@revbalign, $OUT, 1);

exit;
sub PrintAlignments{
	my ($aalign, $balign, $fh, $isrev) = @_;
	if($isrev){
		print "Reverse alignment:\n";
		print $fh "Reverse alignment:\n";
	}else{
		print "Forward alignment:\n";
		print $fh "Forward alignment:\n";
	}
	
	print join(',', @{$aalign}) . "\n";
	print $fh join(',', @{$aalign}) . "\n";
	print join (',', @{$balign}) . "\n";
	print $fh join (',', @{$balign}) . "\n";
	
}
		
sub Alignment{
	my ($swmatrix, $xpos, $ypos, $aarray, $barray, $aalign, $balign) = @_;
		
	my $cont = 1;
	if($swmatrix->[$xpos - 1]->[$ypos - 1] == 0){
		$cont = 0; # This alignment failed completely
	}
	
	my $relhighscore= 0;
	my $relxpos = 0;
	my $relypos = 0;
	
	while($cont){
		#Determine the next value for the matrix high score
		for(my $x = $xpos; $x > 0; --$x) {
			
			if($relhighscore < $swmatrix->[$x-1]->[$ypos-1]) {

				$relhighscore = $swmatrix->[$x-1]->[$ypos-1];
				$relxpos = $x - 1;
				$relypos = $ypos - 1;
			}
		}

		for(my $y = $ypos; $y > 0; --$y) {

			if($relhighscore < $swmatrix->[$xpos-1][$y-1]) {

				$relhighscore = $swmatrix->[$xpos-1][$y-1];
				$relxpos = $xpos -1;
				$relypos = $y - 1;
			}
		}
		
		# perfect diagonal match
		if($relxpos == $xpos - 1 && $relypos == $ypos -1){
			unshift(@{$aalign}, $aarray->[$relxpos -1]);
			unshift(@{$balign}, $barray->[$relypos -1]);
		}else{
			# Gap logic
			if($relxpos == $xpos -1 && $relypos != $ypos - 1){
				# Gap for seq A
				# Add the bases up to the end of the gap for seqB
				for(my $y = $ypos -1; $y > $relypos -1; --$y){
					unshift(@{$balign}, $barray->[$relypos - 1]);
				}
				
				# Add dashes to compensate for seqA
				for(my $x = $ypos -1; $x > $relypos; --$x){
					unshift(@{$aalign}, "-");
				}
			}elsif($relxpos != $xpos - 1 && $relypos == $ypos - 1){
				# Gap for seq B
				for(my $x = $xpos - 1; $x > $relxpos - 1; --$x){
					unshift(@{$aalign}, $aarray->[$relxpos - 1]);
				}
				
				for(my $y = $xpos -1; $y > $relxpos; --$y){
					unshift(@{$balign}, "-");
				}
			}
		
		}	
		$cont = Alignment($swmatrix, $relxpos, $relypos, $aarray, $barray, $aalign, $balign);
	}
	return $cont;
}

sub FindMaxValue{
	my ($swmatrix) = @_;
	
	my $highscore = 0;
	my $xpos = 0;
	my $ypos = 0;
	for(my $x = 0; $x < scalar(@{$swmatrix}); $x++){
		for(my $y = 0; $y < scalar(@{$swmatrix}); $y++){
			if($swmatrix->[$x]->[$y] > $highscore){
				$highscore = $swmatrix->[$x]->[$y];
				$xpos = $x;
				$ypos = $y;
			}
		}
	}
	
	return ($xpos, $ypos, $highscore);
}

sub CreateSWMatrix{
	my ($aarray, $barray) = @_;
	# A is the vertical
	# B is the horizontal
	
	# Initialize the matrix
	my @swmatrix;
	for(my $x = 0; $x < scalar(@{$aarray}); $x++){
		my @tempb;
		for(my $y = 0; $y < scalar(@{$barray}); $y++){
			push(@tempb, 0);
		}
		push(@swmatrix, \@tempb);
	}
	
	# Score the matrix
	for(my $x = 1; $x <= scalar(@{$aarray}); $x++){
		for(my $y = 0; $y <= scalar(@{$barray}); $y++){
			my $compval = 0;
			
			# Position matches
			if($aarray->[$x -1] eq $barray->[$y -1]){
				# Match score = 5
				$compval = $swmatrix[$x-1]->[$y-1] + 5;
			}
			
			$compval = GapPenaltyIncorporation($x, $y, $compval, \@swmatrix);
			
			# Mismatch
			if($aarray->[$x -1] ne $barray->[$y -1]){
				# Mismatch score = -4
				$compval = $swmatrix[$x-1]->[$y-1] - 4;
			}
			
			$compval = GapPenaltyIncorporation($x, $y, $compval, \@swmatrix);
			
			$swmatrix[$x]->[$y] = $compval;
		}
	}
	return @swmatrix;
}

sub GapPenaltyIncorporation{
	my ($x, $y, $compval, $swmatrix) = @_;
	
	# Gap extension penalty = 8
	# No gap open penalty
	for(my $k = $x - 1; $k > 0; --$k){
		if($compval < ($swmatrix->[$x - $k]->[$y] - (8 * $k))){
			$compval = ($swmatrix->[$x - $k]->[$y] - (8 * $k));
		}
	}

	for(my $k = $y - 1; $k > 0; --$k){
		if($compval < ($swmatrix->[$x]->[$y - $k] - (8 * $k))){
			$compval = ($swmatrix->[$x - $k]->[$y] - (8 * $k));
		}
	}

	if($compval < 0){
		$compval = 0;
	}
	return $compval;
}