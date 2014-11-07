#!/usr/bin/perl
# this script is designed to create a table out of merged CNVRs and to cross reference against the copy number indicies
# of animals run through alkan's pipeline. Input is two text files ((1) cnv coords, (2) cn lists) with the same starting animal names
# Output is a tab file with the cnvrs and average copy number for that region.
# cnv coord files should be named beds (simple format, with the animal name)

use strict;
use threads;
use Getopt::Std;

my $threads = 5;
my $bedtools = "";
#my $bedtools = "/mnt/gliu1_usb/dbickhart/BEDTools-Version-2.10.0/bin";
my $usage = "$0 -i \<cnv coords\> -c \<copy number\>\ninputs are text files with the locations of the files\n";

my %opts;
getopt('ic', \%opts);

unless(defined($opts{i}) && defined($opts{c})){
	print $usage;
	exit;
}

my %an_cnv_coords;
my @file_list;
open (IN, "< $opts{i}");
while(my $line = <IN>){
	chomp $line;
	push(@file_list, $line);
	open (CNV, "< $line");
	my @p_seg = split(/\//, $line);
	my @f_seg = split(/\_/, $p_seg[-1]);
	my $an = $f_seg[0];
	$an =~ tr/[a-z]/[A-Z]/;
	while (my $cnv = <CNV>){
		chomp $cnv;
		$cnv =~ s/\r//g;
		my @seg = split(/\t/, $cnv);
		my @type = split(/\_/, $seg[3]);
		push(@{$an_cnv_coords{$an}->{$seg[0]}}, [$seg[1], $seg[2], $type[1]]);
	}
	close CNV;
}
print "Number of cnv files: " . scalar(@file_list) . "\n";
close IN;

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
print "Bedtools merging...\n";
my $cat = join(' ', @file_list);
system("cat $cat > temp");
open(BED, "mergeBed -i temp -nms |") || print "$!\n";
my @merge_cnvrs;
my $wc = 0;
while (my $line = <BED>){
	chomp $line;
	$line =~ s/\r//g;
	my @seg = split(/\t/, $line);
	my @tags = split(/\;/, $seg[3]);
	my $type;
	my $m_string;
	my %unique_tags;
	foreach my $t (@tags){
		my @t_segs = split(/\_/, $t);
		if (!defined($type)){
			$type = $t_segs[1];
		} elsif(length($type) > 1 && $type ne $t_segs[1]){
			$type = "both";
		}
		$t_segs[0] =~ tr/[a-z]/[A-Z]/;
		$unique_tags{$t_segs[0]} = 1;
	}
	$m_string = join(';', sort{$a cmp $b} keys(%unique_tags));
	push (@merge_cnvrs, [$seg[0], $seg[1], $seg[2], $m_string, $type]);
	$wc++;
}
close BED;
print "Merger output lines: $wc\n";

my @animal_list = sort{$a cmp $b} (keys %an_cnv_coords);
my $num = scalar(@animal_list);
my $counter = 1;

open (OUT, "> cnvr_cn_table_buf.tab");
print OUT "cnvr\#\tchr\tstart\tend\t\#animals\tanimals\t";
foreach my $a (@animal_list){
	print OUT "$a CN\t";
}
print OUT "\#overlaps\ttype\toverlapping CNVs\n";
my @out_row;
foreach my $cnvr_ref (@merge_cnvrs){
	my $chr = $cnvr_ref->[0];
	my $start = $cnvr_ref->[1];
	my $end = $cnvr_ref->[2];
	my $ans = $cnvr_ref->[3];
	my $a_count = 1;
	while($ans =~ m/\;/g){$a_count++;}
	push(@out_row, ($counter, $chr, $start, $end, $a_count, $ans));
	#push(@out_row, [$counter, $chr, $start, $end, $a_count, $ans]);

	my @list_of_threads;
	my %cn_results = ();
	my %an_overlap = ();

	foreach my $a_name (@animal_list){
		my ($r_an, $r_val) = cross_index_cn($chr, $start, $end, $an_cn_numbers{$a_name}, $a_name);
		$cn_results{$r_an} = $r_val;
		my ($c_an, $c_use, $c_array) = find_cnv_overlap($chr, $start, $end, $an_cnv_coords{$a_name}, $a_name);
		if ($c_use){
			$an_overlap{$c_an} = $c_array;
		}
	}
	$wc--;
	print " $wc\n";
	foreach my $t_an (@animal_list){
		push(@out_row, $cn_results{$t_an});
	}
	my @temp;
	foreach my $t_an (@animal_list){
		if(exists($an_overlap{$t_an})){
			push(@temp, @{$an_overlap{$t_an}});
		}
	}
	my $over_an = scalar(@temp);
	push(@out_row, $over_an);
	push(@out_row, $cnvr_ref->[4]);
	push(@out_row, @temp);
	
	my $out = join("\t", @out_row);
	print OUT "$out\n";
	@out_row = ();
	$counter++;

}

#system("rm temp");
close OUT;
exit;

sub cross_index_cn (){
	# Takes input coordinates and a reference to a copy number array and calculates average copy number for those coords
	my ($chr, $start, $end, $cn_ref, $a_name) = @_;
	print "$a_name ";
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

sub find_cnv_overlap(){
	# Takes input cnv coordinates and a cnvr coordinate and tries to find overlap in that animal
	my ($chr, $start, $end, $h_ref, $a_name) = @_;
	my $r_use = 0;
	my @r_array;
	foreach my $row (@{$h_ref->{$chr}}){
		if ($row->[1] > $start && $row->[0] < $end){
			push (@r_array, "$a_name\($chr\:$row->[0]\-$row->[1]\)\[$row->[2]\]");
			$r_use = 1;
		}
	}
	return($a_name, $r_use, \@r_array);
}