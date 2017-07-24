#!/usr/bin/perl
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=1
# This is a one-shot script designed to map and order SNP probes from the recombination map to a new assembly
# Output: base.tab <- snpnames and mapping coordinates
# Output: base.segs <- identified chromosome segments
# Output: base.conflicts <- conflicting segments
# Output: base.stats <- general statistics

use strict;
use Getopt::Std;

my $usage = "perl $0 (-a <assembly fasta> -p <probe fasta>) || (-n <nucmer aligns>) -o <output file basename>\n";
my %opts;
getopt('apon', \%opts);

unless(((defined($opts{'a'}) && defined($opts{'p'})) || defined($opts{'n'})) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $unmaps = 0; my $maps = 0; 
my %aligns; # {ochr}->{opos} = [probe, chr, pos, orient]

if(defined($opts{'a'})){
open(my $IN, "module load bwa; bwa mem $opts{a} $opts{p} |") || die "Could not begin BWA alignments!\n";
while(my $line = <$IN>){
	if($line =~ /^@/){
		next;
	}

	chomp $line;
	my @segs = split(/\t/, $line);
	my @rnsegs = split(/\./, $segs[0]);
	if($segs[1] & 2048){
		next;
	}
	if($rnsegs[1] == 0){next;} # Takes care of probes without prior chromosome alignments
	if($segs[2] eq "*"){
		$unmaps++; # Count unmapped probes
	}else{
		$maps++;
	}
	my $orient = ($segs[1] & 16)? "-" : "+";
	$aligns{$rnsegs[1]}->{$rnsegs[2]} = [$rnsegs[0], $segs[2], $segs[3], $orient];
}
close $IN;
}elsif(defined($opts{'n'})){
	open(my $IN, "< $opts{n}") || die "Could not open nucmer aligns!\n";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$maps++;
		my $orient = ($segs[2] > $segs[3])? "-" : "+";
		$aligns{$segs[11]}->{$segs[0]} = ["none", $segs[12], $segs[2], $orient];
	}
	close $IN;
}

open(my $OUT, "> $opts{o}.tab");
open(my $STATS, "> $opts{o}.stats");
open(my $SEGS, "> $opts{o}.segs");
print {$STATS} "Mapping probes: $maps\tUnmapped probes: $unmaps\n";
foreach my $chr (sort{$a <=> $b} keys(%aligns)){
	my ($consensus, $values) = determineConsensus($aligns{$chr});
	print {$STATS} "Ref $chr consensus:";
	for (my $x = 0; $x < scalar(@{$consensus}); $x++){
		print {$STATS} " $consensus->[$x]:$values->[$x]";
	}
	print {$STATS} "\n";

	my ($refblocks, $qblocks) = identifyAndCondenseSegs($aligns{$chr}, $consensus->[0]);
	for(my $x = 0; $x < scalar(@{$refblocks}); $x++){
		my $ref = $refblocks->[$x];
		my $query = $qblocks->[$x];
		my $rlen = abs($ref->[1] - $ref->[0]);
		my $qlen = abs($query->[0] - $query->[1]);
		my $orient = ($query->[0] < $query->[1])? "+" : "-";
		print {$SEGS} "$chr\t$ref->[0]\t$ref->[1]\t$query->[2]\t$query->[0]\t$query->[1]\t$orient\t$rlen\t$qlen\n";
	}

	foreach my $pos (sort{$a <=> $b} keys(%{$aligns{$chr}})){
		my $arrayref = $aligns{$chr}->{$pos};
		print {$OUT} join("\t", @{$arrayref});
		print {$OUT} "\t$chr\t$pos\n";
	}
}
close $OUT;

exit;

sub identifyAndCondenseSegs{
	my ($hashref, $consensus) = @_;
	# Logic: tolerate one deviation in consensus, otherwise condense region into a block
	my @refblock; # [start, end]
	my @queryblock; # [start, end, chr, testbit]

	my @buff; my $count = 0; my $segs = 0; my $skip = 0;
	foreach my $pos (sort{$a <=> $b} keys(%{$hashref})){
		my $query = $hashref->{$pos};
		if($query->[1] eq "*"){next;} # skip unmapped segs
		push(@buff, [$pos, $query]);
		if($count < 2){
			$count++;
			# Fill the initial container buffer
			next;
		}	
		
		if(scalar(@refblock) - 1 < $segs){
			# starting a new segment
			push(@refblock, [$buff[0]->[0], $buff[0]->[0]]);
			push(@queryblock, [$buff[0]->[1]->[2], $buff[0]->[1]->[2], $buff[0]->[1]->[1]]);
		}

		# Test two consecutive probes in the middle of the window to see if they match expectations
		my $comparator = $buff[0]->[1];
		my $test1 = $buff[1]->[1];
		my $test2 = $buff[2]->[1];
		# avg pairwise distance between reference probes in this view
		my $refDist = (($buff[1]->[0] - $buff[0]->[0]) + ($buff[2]->[0] - $buff[1]->[0])) / 2;
		my $t1dist = abs($test1->[2] - $comparator->[2]);
		my $t2dist = abs($test2->[2] - $comparator->[2]);
		
		if($skip){
			$skip = 0;
		}else{
			# passed the test bit for singleton deviations in consensus
			if(($t1dist > 5 * $refDist && $t2dist > 5 * $refDist) ||
				($test1->[1] ne $consensus && $test2->[1] ne $consensus)){
				$segs++; # The conditional now knows to start a new segment
			}elsif($t1dist > 5 * $refDist || $test1->[1] ne $consensus){
				# We don't want singletons to screw up our segments
				$skip = 1;
			}
			# Update the current segments
			$refblock[-1]->[1] = $buff[0]->[0];
               		$queryblock[-1]->[1] = $buff[0]->[1]->[2];
		}

		shift(@buff); # Remove the preceeding buffer item
	}
	if(scalar(@buff) < 3 && scalar(@refblock) < 1){
		# For alignments with fewer lines than chromosomes
		push(@refblock, [$buff[0]->[0], $buff[0]->[0]]);
		push(@queryblock, [$buff[0]->[1]->[2], $buff[0]->[1]->[2], $buff[0]->[1]->[1]]);
		if(scalar(@buff) > 1){
			if($buff[1]->[1]->[2] eq $consensus){
				$refblock[0]->[1] = $buff[1]->[0];
				$queryblock[0]->[1] = $buff[1]->[1]->[2];
			}else{
				push(@refblock, [$buff[1]->[0], $buff[1]->[0]]);
				push(@queryblock, [$buff[1]->[1]->[2], $buff[1]->[1]->[2], $buff[1]->[1]->[1]]);
			}
		}
	}else{
	# Update the final segments for this chr
	$refblock[-1]->[1] = $buff[-1]->[0];
	$queryblock[-1]->[1] = $buff[-1]->[1]->[2];
	}

	return \@refblock, \@queryblock;
}

		

sub determineConsensus{
	my ($hashref) = @_;
	# The input is all of the mapped probes from a reference chr
	# All we need to do is to determine the highest mapping percentile chr from the mapping chr
	my %chrs;
	foreach my $pos (keys(%{$hashref})){
		$chrs{$hashref->{$pos}->[1]} += 1;
	}
	my @consensus = sort{$chrs{$b} <=> $chrs{$a}} keys(%chrs);
	my @values = map{$chrs{$_}} @consensus;
	return \@consensus, \@values;
}
