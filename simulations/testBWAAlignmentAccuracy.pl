#!/usr/bin/perl
# This program is designed to work with my notes on profiling BWA run accuracy

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -i <input sorted bam> -t <input read names> -o <output stats file>\n";

getopt('ito', \%opts);

unless(defined($opts{'i'}) && defined($opts{'t'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# Read name format: EASE:(placeholder):(placeholder):(readnumber)
# Paul's format: (chr)_(position)_?_(read frag number)_?_(read len)_(orient)

# open mapping key
my %chrkey; # key->{read num} -> [read 1, 2]
my %poskey; # key->{read num} -> [read 1, 2] 

my $rnum = 1;
open(my $IN, "< $opts{t}") || die "Could not open key file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	foreach my $s (@segs){
		my @dsegs = split(/_/, $s);
		push(@{$chrkey{$rnum}}, $dsegs[0]);
		push(@{$poskey{$rnum}}, $dsegs[1]);
	}
}
close $IN;

# Now to read the bam and tabulate the scores
my @types = ("CORRECT", "ONECORRECT", "MISSED");
my %readstates; # key->{read num} = value

open(my $IN, "samtools view $opts{i} |");
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my @rsegs = split(/:/, $segs[0]);
	if($segs[1] & 0x800){
		# Split read alignment -- skip!
		next;
	}
	my $fragnum = ($segs[1] & 0x40)? 0 : 1; # 0 for read 1, 1 for read 2
	my $orient = ($segs[1] & 0x10)? "R" : "F";
	
	my $expectC = $chrkey{$rsegs[-1]}->[$fragnum];
	my $expectP = $poskey{$rsegs[-1]}->[$fragnum];
	
	if($expectC ne $segs[2] || $expectP != $segs[3]){
		# shown alignment doesn't match
		if(exists($readstates{$rsegs[-1]})){
			if($readstates{$rsegs[-1]} eq "CORRECT"){
				$readstates{$rsegs[-1]} = "ONECORRECT";
			}
		}else{
			$readstates{$rsegs[-1]} = "MISSED";
		}
	}else{
		if(exists($readstates{$rsegs[-1]})){
			if($readstates{$rsegs[-1]} eq "MISSED"){
				$readstates{$rsegs[-1]} = "ONECORRECT";
			}
		}else{
			$readstates{$rsegs[-1]} = "CORRECT";
		}
	}
}

close $IN;

# Now to do the counting
my %endstats;
my $totalreads = 0;
foreach my $rnum (keys(%readstates)){
	$endstats{$readstates{$rnum}} += 1;
	$totalreads++;
}

# Printing out results
open(my $OUT, "> $opts{o}");
print {$OUT} "Class\tRead Count\tPercentage\n";
foreach my $t (@types){
	if(exists($endstats{$t})){
		my $num = $endstats{$t};
		my $perc = sprintf("%0.4f", ($num / $totalreads));
		print "$t\t$num\t$perc\n";
		print {$OUT} "$t\t$num\t$perc\n";
	}else{
		print "$t\t0\t0.0000\n";
		print {$OUT} "$t\t0\t0.0000\n";
	}
}
close $OUT;
exit;