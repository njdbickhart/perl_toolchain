#!/usr/bin/perl
# This script takes contig names from several assemblies and attempts to align them using a 
# modified Smith-Waterman alignment algorithm. No Affine penalty.

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -a <comparison file one> -b <comparison file two> -o <output alignment>
	Optional: -d	Start the matrix at the lower diagonal rather than at the highest alignment score
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
my $diag = ($opts{'d'})? 1 : 0;

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

my ($forwardmatrix, $ipos, $jpos, $highscore) = CreateSWMatrix(\@adata, \@bdata, $diag);
#my ($xpos, $ypos, $highscore) = FindMaxValue(\@forwardmatrix);
my ($foraalign, $forbalign) = Alignment($forwardmatrix, $ipos, $jpos, \@adata, \@bdata);
open(my $OUT, "> $opts{o}");
print "FileA (top): $opts{a}\n";
print "FileB (bottom): $opts{b}\n";
print $OUT "FileA (top): $opts{a}\n";
print $OUT "FileB (bottom): $opts{b}\n";
PrintAlignments($foraalign, $forbalign, $OUT, $highscore, 0);

my @revb = reverse(@bdata);

my ($reversematrix, $ipos, $jpos, $highscore) = CreateSWMatrix(\@adata, \@revb, $diag);
#($xpos, $ypos, $highscore) = FindMaxValue(\@reversematrix);
my ($revaalign, $revbalign) = Alignment($reversematrix, $ipos, $jpos, \@adata, \@revb);
PrintAlignments($revaalign, $revbalign, $OUT, $highscore, 1);

exit;
sub PrintAlignments{
	my ($aalign, $balign, $fh, $highscore, $isrev) = @_;
	
	if($isrev){
		print "Reverse alignment: $highscore\n";
		print $fh "Reverse alignment: $highscore\n";
	}else{
		print "Forward alignment: $highscore\n";
		print $fh "Forward alignment: $highscore\n";
	}
	
	print join(',', @{$aalign}) . "\n";
	print $fh join(',', @{$aalign}) . "\n";
	print join (',', @{$balign}) . "\n";
	print $fh join (',', @{$balign}) . "\n";
	
}
		
sub Alignment{
	my ($swmatrix, $ipos, $jpos, $aarray, $barray) = @_;
		
	my (@align1, @align2);
	while (1) { 
		if($swmatrix->[$ipos]->[$jpos]->{pointer} eq "none"){
			last;
		}
		
		if ($swmatrix->[$ipos]->[$jpos]->{pointer} eq "diagonal") { 
			push(@align1, $aarray->[$jpos -1]);
			push(@align2, $barray->[$ipos - 1]); 
			$ipos--; 
			$jpos--; 
		} elsif ($swmatrix->[$ipos]->[$jpos]->{pointer} eq "left") { 
			my $len = length($aarray->[$jpos -1]);
			push(@align1, $aarray->[$jpos - 1]);
			
			if($len -2 <= 0){
				push(@align2, "-");
			}else{
				my $filler = "[" . ("-" x ($len - 2)) . "]";
				push(@align2, $filler);
			}
 
			$jpos--; 
		} elsif ($swmatrix->[$ipos]->[$jpos]->{pointer} eq "up") { 
			my $len = length($barray->[$ipos - 1]);
			if($len -1 <= 0){
				push(@align1, "-");
			}else{
				my $filler = "[" . ("-" x ($len - 2)) . "]";
				push(@align1, $filler);
			}
			
			push(@align2, $barray->[$ipos - 1]);
			$ipos--; 
		} 
	}
	@align1 = reverse(@align1);
	@align2 = reverse(@align2);
	
	return \@align1, \@align2;
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
	my ($aarray, $barray, $diag) = @_;
	# A is the vertical
	# B is the horizontal
	my $match = 1;
	my $mismatch = 0;
	my $gap = -1;
	
	# Initialize the matrix
	my @swmatrix;
	$swmatrix[0]->[0]->{score} = 0; 
	$swmatrix[0]->[0]->{pointer} = "none"; 
	for(my $j = 1; $j <= scalar(@{$aarray}); $j++) { 
		$swmatrix[0]->[$j]->{score} = 0; 
		$swmatrix[0]->[$j]->{pointer} = "none"; 
	} 
	for (my $i = 1; $i <= scalar(@{$barray}); $i++) { 
		$swmatrix[$i]->[0]->{score} = 0; 
		$swmatrix[$i]->[0]->{pointer} = "none"; 
	}
	
	# Score the matrix
	my ($ipos, $jpos, $highscore);
	for(my $i = 1; $i <= scalar(@{$barray}); $i++) { 
		for(my $j = 1; $j <= scalar(@{$aarray}); $j++) { 
			my ($diagonal_score, $left_score, $up_score); 
			# calculate match score 
			my $contig1 = $aarray->[$j - 1]; 
			my $contig2 = $barray->[$i - 1]; 
			if ($contig1 eq $contig2) { 
				$diagonal_score = $swmatrix[$i-1]->[$j-1]->{score} + $match; 
			} else { 
				$diagonal_score = $swmatrix[$i-1]->[$j-1]->{score} + $mismatch; 
			} 
			# calculate gap scores 
			$up_score = $swmatrix[$i-1]->[$j]->{score} + $gap; 
			$left_score = $swmatrix[$i]->[$j-1]->{score} + $gap; 
			if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) { 
				$swmatrix[$i]->[$j]->{score} = 0; 
				$swmatrix[$i]->[$j]->{pointer} = "none"; 
				next; 
				# terminate this iteration of the loop 
			}
			# choose best score 
			if ($diagonal_score >= $up_score) { 
				if ($diagonal_score >= $left_score) { 
					$swmatrix[$i]->[$j]->{score} = $diagonal_score; 
					$swmatrix[$i]->[$j]->{pointer} = "diagonal"; 
				} else { 
					$swmatrix[$i]->[$j]->{score} = $left_score; 
					$swmatrix[$i]->[$j]->{pointer} = "left"; 
				} 
			} else { 
				if ($up_score >= $left_score) { 
					$swmatrix[$i]->[$j]->{score} = $up_score; 
					$swmatrix[$i]->[$j]->{pointer} = "up"; 
				} else { 
					$swmatrix[$i]->[$j]->{score} = $left_score; 
					$swmatrix[$i]->[$j]->{pointer} = "left"; 
				} 
			} 
			# set maximum score 
			if ($swmatrix[$i]->[$j]->{score} > $highscore) { 
				$ipos = $i; 
				$jpos = $j; 
				$highscore = $swmatrix[$i]->[$j]->{score}; 
			}
		}
	}
	
	# Useful option if a "more global" alignment is wanted
	if($diag){
		$ipos = scalar(@{$barray}) - 1;
		$jpos = scalar(@{$aarray}) - 1;
	}
	
	return \@swmatrix, $ipos, $jpos, $highscore;
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