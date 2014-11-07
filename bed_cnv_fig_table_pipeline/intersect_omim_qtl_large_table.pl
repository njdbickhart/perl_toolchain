#!/usr/bin/perl
# This script is designed to take the omim.txt and qtl.txt file and generate bed coordinates from it
# it then intersects those bed coordinates with each animal's CNV coordinates and generates a .tab file

use strict;

my $omim = "omim.txt";
my $refseq = "refseq_names.bed";
my $qtl = "qtl.txt";

open(REF, "< $refseq") || die "Could not find $refseq\n";
my %ref_coords;
while(my $line = <REF>){
	chomp $line;
	my @segs = split(/\t/, $line);
	push (@{$ref_coords{$segs[3]}}, ($segs[0], $segs[1], $segs[2]));
}
close REF;
open (QTL, " < $qtl") || die "Could not find $qtl\n";
my %qtl_sorter;
my @qtl_data;
while(my $line = <QTL>){
	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line);
	my $c = shift(@segs);
	push(@{$qtl_sorter{$c}}, [@segs]);
}
foreach my $c (sort {$a cmp $b} keys(%qtl_sorter)){
	foreach my $a (@{$qtl_sorter{$c}}){
		push (@qtl_data, [($c, $a->[0], $a->[1], $a->[2], $a->[3])]);
	}
}
close QTL;

open (OMIM, "< $omim") || die "Could not find $omim\n";
my $header = <OMIM>;
my %omim_coords;
while(my $line = <OMIM>){
	chomp $line;
	my @temp_coords = ();
	my @segs = split(/\t/, $line);
	my @ref = split(/\,\s/, $segs[1]);
	my $v_coords;
	foreach my $r (@ref){
		if(exists($ref_coords{$r})){
			$v_coords = $ref_coords{$r};
		}
		if (scalar(@temp_coords > 0) && exists($ref_coords{$r})){
			$temp_coords[2] = $v_coords->[2];
		}elsif(scalar(@temp_coords == 0) && exists($ref_coords{$r})){
			push(@temp_coords, @{$v_coords});
		}
	}
	if (scalar(@temp_coords) > 0){
		push(@temp_coords, $segs[0], $segs[1]);
		my $chr = shift(@temp_coords);
		push(@{$omim_coords{$chr}}, [@temp_coords]);
	}
}
close OMIM;
open (OUT, "> omim_qlt_comb_table.tab");
print OUT "No\tchr\tstart\tend\tlength\trefseq symbols\tOMIM No\tOMIM information\tQTL No\tQTL information\n";
my $count = 1;
foreach my $q (@qtl_data){
	my $chr = $q->[0];
	my $start = $q->[1];
	my $end = $q->[2];
	my $len = $end - $start;
	my $rseq = "";
	my $omim_num = 0;
	my $omim_inf = "";
	my $qtl_count = $q->[3];
	my $qtl_inf = $q->[4];
	my @temp_vals = ("", ""); #0: annotation, 1: refseq
	foreach my $om (@{$omim_coords{$chr}}){
		if ($om->[0] <= $end && $om->[1] >= $start){
			$temp_vals[0] .= "$om->[2];";
			$temp_vals[1] .= "$om->[3], ";
		}
	}
	if (length($temp_vals[0]) > 1){
		chop $temp_vals[0];
		chop $temp_vals[1]; chop $temp_vals[1];
		my @om_array = split(/\;/, $temp_vals[0]);
		$omim_num = scalar(@om_array);
		$rseq = $temp_vals[1];
		$omim_inf = $temp_vals[0];
	}
	print OUT "$count\t$chr\t$start\t$end\t$len\t$rseq\t$omim_num\t$omim_inf\t$qtl_count\t$qtl_inf\n";
	$count++;
}
close OUT;
exit;