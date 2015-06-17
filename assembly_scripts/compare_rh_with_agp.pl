#!/usr/bin/perl
# This script is designed to pull problem contigs from both the RH and AGP maps to try to get quick
# location comparisons for conflict resolution

use strict;
use Getopt::Std;

my $usage = "perl $0 -a <agp file> -r <rh file> -c <chromosome> -p <contig> -l <lines to pull above and below> -s <skip instance count [default: none]>\n";
my %opts;
getopt('arcpls', \%opts);
# OK, the goal is to condense the RH map down to contitutive parts

unless(defined($opts{'a'}) && defined($opts{'r'})){
	print $usage;
	exit;
}

my @agpdata; # [lines]-> [name, start, end, orient, supscaff]
my @rhdata; # [lines]-> [name, rhstart, rhend, start, end]
# Start with agp file first
my $agpbuffer = $opts{'l'};
my @linebuff;
my $skipthresh = (defined($opts{'s'}))? $opts{'s'} : 0;
my $skipcount = 0;
my $agp;
open($agp, "< $opts{a}") || die "Could not open $opts{a} agp file!\n";
while(my $line = <$agp>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if($segs[5] eq $opts{'p'} && $segs[0] eq $opts{'c'} && $skipcount == $skipthresh){
		# Found the first instance! Working with it
		# starting with previous data
		foreach my $r (@linebuff){
			push(@agpdata, [$r->[5], $r->[6], $r->[7], $r->[8], $r->[9]]);
		}
		push(@agpdata, ["\[$segs[5]\]", $segs[6], $segs[7], $segs[8], $segs[9]]);
		for(my $x = 0; $x < $agpbuffer; $x++){
			my $l = <$agp>;
			chomp $l;
			@segs = split(/\t/, $l);
			if($segs[5] eq $opts{'p'}){
				$segs[5] = "\[$segs[5]\]";
			}
			push(@agpdata, [$segs[5], $segs[6], $segs[7], $segs[8], $segs[9]]);
		}
		last; # We're done
	}elsif($segs[5] eq $opts{'p'} && $segs[0] eq $opts{'c'} && $skipcount != $skipthresh){
		$skipcount++;
	}
	if(scalar(@linebuff) >= $agpbuffer){
		shift(@linebuff);
	}
	push(@linebuff, \@segs);
}
close $agp;
@linebuff = (); # clearing the buffer
$skipcount = 0;

# Now for the RH file -- this will be tricky because of multiple listings
my $rhbuff = $opts{'l'} * 60;
#my @rhdata; # [lines]-> [name, rhstart, rhend, start, end]
open(my $rh, "< $opts{r}") || die "Could not open $opts{r} rh file\n";
while(my $line = <$rh>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my @qsegs = split(/\_/, $segs[5]);
	$segs[5] = $qsegs[0];
	if($segs[1] eq $opts{'c'} && $qsegs[0] eq $opts{'p'} && $skipcount == $skipthresh){
		# Found the start! Now to progress back through the buffer to get the lines we need to pull
		my @order;
		foreach my $r (@linebuff){
			if(scalar(@order) == 0){
				push(@order, $r->[5]);
			}
			if($r->[5] ne $order[-1]){
				push(@order, $r->[5]);
			}
		}
		
		# OK, let's count how many we need to grab
		my $startorder = 0;
		if(scalar(@order) > $opts{'l'}){
			$startorder = scalar(@order) - $opts{'l'}; 
		}
		my $startcontig = $order[$startorder];
		
		my $begin = 0;
		for(my $x = 0; $x < scalar(@linebuff); $x++){
			if($begin || $linebuff[$x]->[5] eq $startcontig){
				if(scalar(@rhdata) == 0){
					push(@rhdata, [$linebuff[$x]->[5], $linebuff[$x]->[2], 0, $linebuff[$x]->[6], 0]);
					$begin = 1;
					next;
				}
				if($rhdata[-1]->[0] eq $linebuff[$x]->[5]){
					$rhdata[-1]->[2] = $linebuff[$x]->[2];
					$rhdata[-1]->[4] = $linebuff[$x]->[6];
				}
				else{
					push(@rhdata, [$linebuff[$x]->[5], $linebuff[$x]->[2], 0, $linebuff[$x]->[6], 0]);
				}
				$begin = 1;
			}
		}
		$segs[5] = "\[$segs[5]\]";
		push(@rhdata, [$segs[5], $segs[2], 0, $segs[6], 0]);
		
		# Now to feed forward to get the rest of the lines
		my $lcount = 0;
		while(my $l = <$rh>){
			chomp $l;
			my @segs = split(/\t/, $l);
			@qsegs = split(/\_/, $segs[5]);
			$segs[5] = $qsegs[0];
			if($segs[5] eq $opts{'p'}){
				$segs[5] = "\[$segs[5]\]";
			}
			if($rhdata[-1]->[0] eq $segs[5]){
				$rhdata[-1]->[2] = $segs[2];
				$rhdata[-1]->[4] = $segs[6];
			}
			else{
				push(@rhdata, [$segs[5], $segs[2], 0, $segs[6], 0]);
				$lcount++;
			}
			if($lcount >= $opts{'l'}){
				last;
			}
		}
		last;
	}elsif($segs[1] eq $opts{'c'} && $qsegs[0] eq $opts{'p'} && $skipcount != $skipthresh){
		# Accounts for single markers being missed -- also prevents long lines from saturating
		# the skipcounter
		if($linebuff[-1]->[5] eq $qsegs[0]){
			
		}else{
			$skipcount++;
		}
	}
			
	
	if(scalar(@linebuff) >= $rhbuff){
		shift(@linebuff);
	}
	push(@linebuff, \@segs);
}

close $rh;

# OK, I know that the two data arrays might be different sizes, so let's try to plan this out
my $agpnum = scalar(@agpdata);
my $rhnum = scalar(@rhdata);
my $largest = ($agpnum > $rhnum)? $agpnum : $rhnum;

#my @agpdata; # [lines]-> [name, start, end, orient, supscaff]
#my @rhdata; # [lines]-> [name, rhstart, rhend, start, end]
print "AGPcontig\tAGPcRange\tAGPOrient\tAGPsupscaff\tRHcontig\tRHrange\tRHcRange\n";
for(my $x = 0; $x < $largest; $x++){
	if($x < scalar(@agpdata)){
		print $agpdata[$x]->[0] . "\t" . $agpdata[$x]->[1] . "-" . $agpdata[$x]->[2] . "\t";
		print $agpdata[$x]->[3] . "\t" . $agpdata[$x]->[4] . "\t";
	}else{
		print "\t\t\t\t";
	}
	
	if($x < scalar(@rhdata)){
		# Account for single data point alignments
		my $rrange = ($rhdata[$x]->[2] == 0)? $rhdata[$x]->[1] : $rhdata[$x]->[1] . "-" . $rhdata[$x]->[2];
		my $crange = ($rhdata[$x]->[4] == 0)? $rhdata[$x]->[3] : $rhdata[$x]->[3] . "-" . $rhdata[$x]->[4];
		
		print $rhdata[$x]->[0] . "\t$rrange\t$crange\n";
	}else{
		print "\t\t\n";
	}
}

exit;