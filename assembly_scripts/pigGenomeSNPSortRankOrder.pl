#!/usr/bin/perl
# This script ranks SNP probes from the different assemblies in order of their supposed position on the chromosomes
# Alt haplotypes and unmapped probes are precluded from analysis by ranking as -1

use strict;
my $usage = "perl $0 <input tab files> <output ranking file>\n";

unless(scalar(@ARGV) == 2){
        print $usage;
        exit;
}

chomp(@ARGV);

my %ss10 = (
"NC_010461.4" => 19,
"NC_010462.2" => 20,
"NC_010443.4" => 1,
"NC_010444.3" => 2,
"NC_010445.3" => 3,
"NC_010446.4" => 4,
"NC_010447.4" => 5,
"NC_010448.3" => 6,
"NC_010449.4" => 7,
"NC_010450.3" => 8,
"NC_010451.3" => 9,
"NC_010452.3" => 10,
"NC_010453.4" => 11,
"NC_010454.3" => 12,
"NC_010455.4" => 13,
"NC_010456.4" => 14,
"NC_010457.4" => 15,
"NC_010458.3" => 16,
"NC_010459.4" => 17,
"NC_010460.3" => 18, 
"other" => 21);

my %ros = (
"NC_010461.5" => 19,
"NC_010462.3" => 20,
"NC_010443.5" => 1,
"NC_010444.4" => 2,
"NC_010445.4" => 3,
"NC_010446.5" => 4,
"NC_010447.5" => 5,
"NC_010448.4" => 6,
"NC_010449.5" => 7,
"NC_010450.4" => 8,
"NC_010451.4" => 9,
"NC_010452.4" => 10,
"NC_010453.5" => 11,
"NC_010454.4" => 12,
"NC_010455.5" => 13,
"NC_010456.5" => 14,
"NC_010457.5" => 15,
"NC_010458.4" => 16,
"NC_010459.5" => 17,
"NC_010460.4" => 18,
"other" => 21);

my %marc = (
"CM009104.1" => 19,
"NPJO01000144.1" => 20,
"CM009086.1" => 1,
"CM009087.1" => 2,
"CM009088.1" => 3,
"CM009089.1" => 4,
"CM009090.1" => 5,
"CM009091.1" => 6,
"CM009092.1" => 7,
"CM009093.1" => 8,
"CM009094.1" => 9,
"CM009095.1" => 10,
"CM009096.1" => 11,
"CM009097.1" => 12,
"CM009098.1" => 13,
"CM009099.1" => 14,
"CM009100.1" => 15,
"CM009101.1" => 16,
"CM009102.1" => 17,
"CM009103.1" => 18, 
"other" => 21);

# OK, so the file columns should be: 0. probe, 1. refchr, 2. refpos, 3. rchr, 4. rpos, 5. schr, 6. spos, 7. mchr, 8. mpos
my %sorter = (
        "ROS" => \%ros,
        "SS10" => \%ss10,
        "MARC" => \%marc);
my %tabfields = (
	"ROS" => 3,
	"SS10" => 5,
	"MARC" => 7);
my %probes; # {probe} -> [ref rank, rrank, srank, mrank]
my %datasets; #{data} -> [] ->[probe, chr, pos]

open(my $IN, "< $ARGV[0]");
my $counter = 1;
while(my $line = <$IN>){
        chomp $line;
        my @segs = split(/\t/, $line);
        $probes{$segs[0]} = [$counter];
        $counter++;
        foreach my $id (qw(ROS SS10 MARC)){
        	my $num = $tabfields{$id};
        	my %key = %{$sorter{$id}};
        	if(!exists($key{$segs[$num]})){
        		# Catch unplaced chromosomes and put them in the "other" category for sorting
        		$segs[$num] = "other";
        	}
        	push(@{$datasets{$id}}, [$segs[0], $segs[$num], $segs[$num + 1]]);
        }
        
        #push(@{$datasets{"ROS"}}, [$segs[0], $segs[3], $segs[4]]);
        #push(@{$datasets{"SS10"}}, [$segs[0], $segs[5], $segs[6]]);
        #push(@{$datasets{"MARC"}}, [$segs[0], $segs[7], $segs[8]]);
}
close $IN;

# Now to sort each one and rerank
my $icol = 1;
foreach my $d (qw(ROS SS10 MARC)){
        my @array = @{$datasets{$d}};
        my %key = %{$sorter{$d}};
        my @sort = sort{$key{$a->[1]} <=> $key{$b->[1]} || $a->[2] <=> $b->[2]} @array;
        my $c = 1;
        foreach my $k (@sort){
                if($k->[1] eq "*"){
                        $probes{$k->[0]}->[$icol] = -1;
                        next;
                }
                $probes{$k->[0]}->[$icol] = $c;
                $c++;
        }
        $icol++;
}

open(my $OUT, "> $ARGV[1]");
foreach my $p (sort{$probes{$a}->[0] <=> $probes{$b}->[0]}keys(%probes)){
        print {$OUT} "$p\t" . join("\t", @{$probes{$p}}) . "\n";
}
close $OUT;

exit;
