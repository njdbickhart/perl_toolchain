#!/usr/bin/perl
# This script takes STDIN and counts bed length intervals to give back several stats
# I got tired of writing the same one-liner all the time, basically
# 12/8/2015:	Migrated to my perl repository -- added stdin and ARGV input

use strict;
my $FH;
if(scalar(@ARGV) == 0){
	$FH = \*STDIN;
}else{
	chomp(@ARGV);
	open($FH, "< $ARGV[0]");
}

my @int_lens;
while (my $line = <STDIN>){
	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line);
	if ($segs[1] > $segs[2]){
		($segs[1], $segs[2]) = ($segs[2], $segs[1]);
	}
	push(@int_lens, ($segs[2] - $segs[1]));
}

my $int_num = scalar(@int_lens);
my $int_tot = total_length(\@int_lens);
my $int_avg = average(\@int_lens);
my $int_stdev = standard_deviation(\@int_lens);
my ($small, $large, $int_median) = median(\@int_lens);
print "\tInterval Numbers:\t$int_num\n";
print "\tTotal Length:\t\t$int_tot\n";
print "\tLength Average:\t\t$int_avg\n";
print "\tLength Median:\t\t$int_median\n";
print "\tLength Stdev:\t\t$int_stdev\n";
print "\tSmallest Length:\t$small\n";
print "\tLargest Length:\t\t$large\n";

exit;

sub standard_deviation {
	my $array_ref = shift(@_);
	my @numbers = @{$array_ref};
	return undef unless(scalar(@numbers));
	
	my $total1 = 0;
	foreach my $num (@numbers) {
		$total1 += $num;
	}
	my $mean1 = $total1 / (scalar @numbers);
	my $total2 = 0;
	foreach my $num (@numbers) {
		$total2 += ($mean1-$num)**2;
	}
	my $mean2 = $total2 / (scalar @numbers);
	my $std_dev = sqrt($mean2);
	return $std_dev;
}
sub average {
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	if (scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}
sub median {
	my $rpole = shift @_;
	my @tpole = @$rpole;
	if (scalar(@tpole) == 0){
	return 0;
	}
	my $ret;

	my @pole = sort {$a <=> $b} @tpole;

	if( (@pole % 2) == 1 ) {
		$ret = $pole[((@pole+1) / 2)-1];
	} else {
		$ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
	}

	return $pole[0], $pole[-1], $ret;
}
sub total_length {
	my $a_ref = shift(@_);
	if (scalar(@{$a_ref}) == 0){
		return 0;
	}
	my $tot = 0;
	foreach my $len (@{$a_ref}){
		$tot += $len;
	}
	return $tot;
}