#!/usr/bin/perl
# This is a one-shot script for checking scaffolding orientations and decisions against the RH map
# Input is an AGP file and the RH map ordering of contigs
# RH map ordering file:
# chr rhpos contig orientation

use strict;

my $usage = "perl $0 <input rh map> <input agp file>\n";

chomp(@ARGV);
unless(scalar(@ARGV) == 2){
	print $usage;
	exit;
}

# read in RH map order assignments
my %probeLookup; # lookup for probe index values
my @probes; # []->[chr rhpos contig orientation]

open(my $IN, "< $ARGV[0]") || die "Could not open RH map positions\n";
my $rhcount = 0;
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	push(@probes, \@segs);
	$probeLookup{$segs[2]} = $rhcount;
	$rhcount++;
}

close $IN;

# read in and store scaffolds
my %scaffolds; # {scaffold num} -> [] ->[contig, orientation]
open(my $IN, "< $ARGV[1]") || die "Could not open agp file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[4] ne "W"){next;}
	$segs[5] =~ s/(utg\d{1,8})_.+/$1/;
	push(@{$scaffolds{$segs[0]}}, [$segs[5], $segs[8]]);
}

close $IN;

my $wrongChr = 0;	# Chromosome mistake
my $notFoundScaff = 0;  # Scaffold was too chimeric to determine orientation
my $orientMistake = 0;	# Scaffold in wrong order or orientation
open(my $LOG, "> $ARGV[1].log");
# Main comparison routine
foreach my $scaff (sort {$a <=> $b} keys(%scaffolds)){
	# Determine consensus chromosome
	my $chr = determineConsensus($scaffolds{$scaff}, \@probes, \%probeLookup);
	
	# Determine the optimal orientation on the target chr
	my ($orient, $found) = determineOptimumOrient($scaffolds{$scaff}, \@probes, \%probeLookup, $chr);
	
	if(!$found){
		print {$LOG} "Short_Chimeric\t$scaff\tnum_components:" . scalar(@{$scaffolds{$scaff}}) . "\n";
		$notFoundScaff++;
		next;
	}
	
	# Now, loop through and count based on expected orientation and chr
	my $lastIdx = -1;
	for(my $x = 0; $x < scalar(@{$scaffolds{$scaff}}); $x++){
		my $row = $scaffolds{$scaff}->[$x];
		if(!exists($probeLookup{$row->[0]})){
			next; # This contig does not exist on the RH map, so we don't count it
		}
		
		my $idx = $probeLookup{$row->[0]};
		my $curRow = $probes[$idx];
		if($lastIdx == -1){
			if($curRow->[0] ne $chr){
				$wrongChr++;
				print {$LOG} "Wrong_chr\t$scaff\t$row->[0]\tconsensus_scaff_chr:$chr\trh_chr:" . $curRow->[0] . "\n";
				next;
			}
			$lastIdx = $idx;
			next;
		}else{			
			my $lastRow = $probes[$lastIdx];
			if($curRow->[0] ne $chr){
				$wrongChr++;
				print {$LOG} "Wrong_chr\t$scaff\t$row->[0]\tconsensus_scaff_chr:$chr\trh_chr:" . $curRow->[0] . "\n";
				$lastIdx = -1; # reset the idx counter to try again
				next;
			}
			
			# probes []->[chr rhpos contig orientation]
			my $sameOrient = ($orient eq "+")? ($curRow->[3] eq $row->[1])? 1 : 0 : ($curRow->[3] ne $row->[1])? 1 : 0;
			# Easy comparison
			unless($lastRow->[2] eq $scaffolds{$scaff}->[$x - 1]->[0]
				&& $sameOrient){
				$orientMistake++;
				print {$LOG} "Orientation_error\t$scaff\t$row->[0]\trh_align_orient:$orient\trh_orient:". $curRow->[3] . "\tscaff_orient:" . $row->[1] . "\n";
			}
		}
		$lastIdx = $idx;
	}
}
	
print "Wrong chr: $wrongChr\n";
print "Orient mistake: $orientMistake\n";
print "Scaffold chimeric: $notFoundScaff\n";
close $LOG;
exit;

sub determineOptimumOrient{
	my ($sref, $pref, $plookup, $chr) = @_;
	
	# Find the first entry that matches the target chr
	my %orients; # To determine best orientation

	my $found = 0;
	for(my $x = 1; $x < scalar(@{$sref}) - 1; $x++){
		my $row = $sref->[$x];
		if(!exists($plookup->{$row->[0]})){
			next;
		}
		
		my $idx = $plookup->{$row->[0]};
		if($pref->[$idx]->[0] eq $chr){
			# Check to see if the next entry is in order
			# Change: removed the orientation check logic to pair as many scaffolds as possible
			my $for = $sref->[$x + 1];
			if(($for->[0] eq $pref->[$idx + 1]->[2])){
				$orients{"+"} += 1;
				$found = 1;
			}elsif($for->[0] eq $pref->[$idx - 1]->[2]){
				$orients{"-"} += 1;
				$found = 1;
			}
		}
	}
	
	my $maxv = 0; my $bestorient = "";
	foreach my $key (keys(%orients)){
		if($orients{$key} > $maxv){
			$maxv = $orients{$key};
			$bestorient = $key;
		}
	}
	
	return($bestorient, $found);
}

sub determineConsensus{
	my ($sref, $pref, $plookup) = @_;
	my %chrs;
	
	foreach my $row (@{$sref}){
		if(!exists($plookup->{$row->[0]})){
			next;
		}
		my $idx = $plookup->{$row->[0]};
		my $chr = $pref->[$idx]->[0];
		$chrs{$chr} += 1;
	}
	
	my $maxv = 0; my $maxchr = "";
	foreach my $key (keys(%chrs)){
		if($chrs{$key} > $maxv){
			$maxv = $chrs{$key};
			$maxchr = $key;
		}
	}
	
	return $maxchr;
}
	