#!/usr/bin/perl
# This script takes individual animal cnv files and calculates the average copynumber for each of them

use strict;
use Getopt::Std;
use Spreadsheet::WriteExcel;

my $bedtools = "/home/derek/share/BEDTools-Version-2.10.1/bin";

my $usage = "$0 -i \<cnv file lists\> -c \<1kb copy number\>
\ninputs are text files with the locations of the files\n";

my %opts;
getopt('icvo', \%opts);

unless(defined($opts{i}) && defined($opts{c})){
	print $usage;
	exit;
}
if (defined($opts{o})){
	print "Output prefix: $opts{o}\n";
} else {
	$opts{o} = "";
}
print "Working on 1kb cn windows\n";
my %an_cn_numbers;
open (IN, "< $opts{c}");
while (my $line = <IN>){
	chomp $line;
	open (CN, "< $line");
	my @p_seg = split(/\//, $line);
	my @f_seg = split(/[._]/, $p_seg[-1]);
	my $an = $f_seg[0];
	$an =~ tr/[a-z]/[A-Z]/;
	while (my $cn = <CN>){
		chomp $cn;
		my @seg = split(/\t/, $cn);
		push(@{$an_cn_numbers{$an}->{$seg[0]}}, [$seg[1], $seg[2], $seg[3]]);
	}
	close CN;
}
close IN;
=pod
print "Working on 5kb cn windows\n";
my %k5_cn_numbers;
open (IN, "< $opts{v}");
while (my $line = <IN>){
	chomp $line;
	open (CN, "< $line");
	my @p_seg = split(/\//, $line);
	my @f_seg = split(/\_/, $p_seg[-1]);
	my $an = $f_seg[0];
	$an =~ tr/[a-z]/[A-Z]/;
	while (my $cn = <CN>){
		chomp $cn;
		my @seg = split(/\t/, $cn);
		push(@{$k5_cn_numbers{$an}->{$seg[0]}}, [$seg[1], $seg[2], $seg[3]]);
	}
	close CN;
}
close IN;
=cut
my %an_cnv_coords;
my @cnv_f_list;
open (IN, "< $opts{i}");
while(my $line = <IN>){
	chomp $line;
	$line =~ s/\r//g;
	push(@cnv_f_list, $line);
	open (CNV, "< $line");
	my @p_seg = split(/\//, $line);
	my @f_seg = split(/\_/, $p_seg[-1]);
	my $an = $f_seg[0];
	$an =~ tr/[a-z]/[A-Z]/;
	while (my $cnv = <CNV>){
		chomp $cnv;
		$cnv =~ s/\r//g;
		my @seg = split(/\t/, $cnv);
		my @event = split(/\_/, $seg[3]);
		push(@{$an_cnv_coords{$an}->{$seg[0]}}, [$seg[1], $seg[2], $event[1]]);
	}
	close CNV;

}
print "Number of cnv files: " . scalar(@cnv_f_list) . "\n";
close IN;


my @animal_list = sort{$a cmp $b} (keys %an_cnv_coords);
my $num = scalar(@animal_list);

foreach my $animal (@animal_list){
	open(OUT, "> $animal$opts{o}\_cnv_cn.bed");
	my $workbook = Spreadsheet::WriteExcel->new("$animal$opts{o}\_cnv_copynumber.xls");
	my $worksheet = $workbook->add_worksheet($animal);
	$worksheet->write(0,0,"chr");
	$worksheet->write(0,1,"start");
	$worksheet->write(0,2,"end");
	$worksheet->write(0,3,"type");
	$worksheet->write(0,4,"1kb avg cn");
	#$worksheet->write(0,5,"5kb avg cn");
	my $row = 1;
	print "Working on $animal\n";
	foreach my $chr (sort {$a cmp $b} keys(%{$an_cnv_coords{$animal}})){
		foreach my $an_row (@{$an_cnv_coords{$animal}->{$chr}}){
			my $s = $an_row->[0];
			my $e = $an_row->[1];
			my $t = $an_row->[2];
			my ($holder, $avg_cn) = cross_index_cn($chr, $s, $e, $an_cn_numbers{$animal}, $animal);
			#my ($holder, $avg_5k) = cross_index_cn($chr, $s, $e, $k5_cn_numbers{$animal}, $animal);
			$worksheet->write($row,0,"$chr");
			$worksheet->write($row,1,"$s");
			$worksheet->write($row,2,"$e");
			$worksheet->write($row,3,"$t");
			$worksheet->write($row,4,"$avg_cn");
			#$worksheet->write($row,5,"$avg_5k");
			print OUT "$chr\t$s\t$e\t$avg_cn\n";
			$row++;
		}
	}
	$workbook->close();
	close OUT;
}


exit;

sub cross_index_cn (){
	# Takes input coordinates and a reference to a copy number array and calculates average copy number for those coords
	my ($chr, $start, $end, $cn_ref, $a_name) = @_;
	#print "$a_name ";
	my @values;
	my $term = 0;
	
	foreach my $l_ref (@{$cn_ref->{$chr}}){
		if ($l_ref->[1] > $start && $l_ref->[0] < $end){
			push (@values, $l_ref->[2]);
			$term = 1;
		}
	}
	
	my $avg1;
	my $avg2;
	my $num = scalar(@values);
	if ($num == 0){
		$num = 1;
	}
	foreach my $v (@values){
		$avg1 += $v;
	}
	$avg2 = $avg1 / $num;
	return($a_name, $avg2);

}