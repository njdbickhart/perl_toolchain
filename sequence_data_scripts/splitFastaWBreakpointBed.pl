#!/usr/bin/perl
# this is a one-shot script designed to create split fasta entries from a bed file with breakpoint coordinates
# 8/25/2015: 	update to account for breakpoints too close to the beginning and too close to the end of the chromosome
#		also updated program to print out samtools faidx coordinates
# 8/26/2015:	Added a minimum size filter criteria to the script. Also added an N ratio filter for split fragments

use strict;
use Getopt::Std;

my $usage = "perl $0 -f <indexed fasta file> -o <output fasta file> -b <bed file with breakpoints> -s <samtools split coordinates> -m <minimum fasta size> -n <ratio of N's to remove from split segment[float]>\n";
my %opts;
getopt('fobsmn', \%opts);

unless(defined($opts{'f'})){
	print $usage;
	exit;
}

my $minsize = 100;
if(defined($opts{'m'})){
	$minsize = $opts{'m'};
}

my $nratio = 1;
if(defined($opts{'n'})){
	$nratio = $opts{'n'};
	if($nratio > 1){
		print STDERR "Error with nratio value! Must be below 1!\n";
		print STDERR $usage;
		exit;
	}
}

unless( -s "$opts{f}.fai"){
	print STDERR "Error! Please index fasta!\n";
	exit;
}

my %coords; # {ctg} -> [] ->[start, end]
my %chrends; # {ctg} = endcoord
# store the bed coordinates
open(my $IN, "< $opts{b}") || die "Could not open bed file!\n";
while(my $line = <$IN>){
	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line);
	push(@{$coords{$segs[0]}}, [$segs[1], $segs[2]]);
}
close $IN;

# Get chr endpoints
open(my $IN, "< $opts{f}.fai") || die "Could not open fasta fai file!\n";
while(my $line = <$IN>){
	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line);
	if(exists($coords{$segs[0]})){
		$chrends{$segs[0]} = $segs[1];
	}
}
close $IN;

open(my $OUT, "> $opts{o}");
open(my $SPLIT, "> $opts{s}");
# Work on split contigs first
foreach  my $ctg (keys(%coords)){
	my @breaks = @{$coords{$ctg}};
	
	if(scalar(@breaks) == 1){
		my $ctr = 1;
		my $start = $breaks[0]->[0];
		my $end = $breaks[0]->[1];
		
		if($start >= 100 && $start - 1 > $minsize){
			if(nCheck($nratio, $opts{'f'}, "$ctg\:1-$start", $start - 1)){
				open(my $IN, "samtools faidx $opts{f} $ctg\:1-$start |");
				my $h = <$IN>;
				print $OUT ">$ctg.$ctr\n";
				print $SPLIT "$ctg.$ctr\t1\t$start\n";
				while(my $line = <$IN>){
					print $OUT $line;
				}
				close $IN;
				$ctr++;
			}
		}
		
		my $chrend = $chrends{$ctg};
		if($end <= $chrend - 100 && $chrend - $end > $minsize){	
			if(nCheck($nratio, $opts{'f'}, "$ctg\:$end-$chrend", $chrend - $end)){
				open(my $IN, "samtools faidx $opts{f} $ctg\:$end-$chrend |");
				my $h = <$IN>;
				print $OUT ">$ctg.$ctr\n";
				print $SPLIT "$ctg.$ctr\t$end\t$chrend\n";
				while(my $line = <$IN>){
					print $OUT $line;
				}
				close $IN;
			}
		}
	}else{
		my $prev = 1;
		my $ctr = 0;
		for(my $x = 0; $x < scalar(@breaks); $x++){
			my $end = $breaks[$x]->[0];
			
			if($end >= 100 && $end - $prev > $minsize){
				if(nCheck($nratio, $opts{'f'}, "$ctg\:$prev-$end", $end - $prev)){
					$ctr++;
					open(my $IN, "samtools faidx $opts{f} $ctg\:$prev-$end |");
					my $h = <$IN>;
					print $OUT ">$ctg.$ctr\n";
					print $SPLIT "$ctg.$ctr\t$prev\t$end\n";
					while(my $line = <$IN>){
						print $OUT $line;
					}
					close $IN;
				}
			}
			
			$prev = $breaks[$x]->[1];
		}
		
		my $chrend = $chrends{$ctg};
		if($prev <= $chrend - 100 && $chrend - $prev > $minsize){
			if(nCheck($nratio, $opts{'f'}, "$ctg\:$prev-$chrend", $chrend - $prev)){
				$ctr++;
				open(my $IN, "samtools faidx $opts{f} $ctg\:$prev-$chrend |");
				my $h = <$IN>;
				print $OUT ">$ctg.$ctr\n";
				print $SPLIT "$ctg.$ctr\t$prev\t$chrend\n";
				while(my $line = <$IN>){
					print $OUT $line;
				}
				close $IN;
			}
		}
	}
}

print STDERR "Working on non-altered fasta headers\n";
# Now, open up the whole fasta file and only print it out if the contig didn't need to be split
open(my $IN, "< $opts{f}");
my $print = 0;
while(my $line = <$IN>){
	chomp $line;
	if($line =~ /^>/){
		my ($ctg) = $line =~ />(.+)/;
		if(!exists($coords{$ctg})){
			$print = 1;
			print $OUT "$line\n";
		}else{
			$print = 0;
		}
	}elsif($print){
		print $OUT "$line\n";
	}
}
close $IN;
close $OUT;
close $SPLIT;

exit;

# Returns "1" if the check passes; "0" if the check fails and the segment should be filtered
sub nCheck{
	my ($nratio, $fasta, $samseg, $seqlen) = @_;
	if($nratio == 1){
		# ncheck is disabled
		return 1;
	}
	my $nsum = 0;
	
	open(my $IN, "samtools faidx $fasta $samseg |");
	my $h = <$IN>; 
	while(my $line = <$IN>){
		my $ns = ($line =~ tr/N/N/);
		$nsum += $ns;
	}
	close $IN;
	
	if($nsum / $seqlen >= $nratio){
		return 0;
	}else{
		return 1;
	}
}