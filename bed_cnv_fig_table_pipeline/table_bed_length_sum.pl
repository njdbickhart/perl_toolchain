#!/usr/bin/perl
# This script takes multiple bed files in @ARGV input and creates a formatted table for stdout
# 12/8/2015:	Migrated to my perl_toolchain repository -- added markdown inspired table

use strict;

my $usage = "$0 [-c <markdown format>] <bed file1> ...\n";

if(scalar(@ARGV) == 0){
	print $usage;
	exit;
}
chomp(@ARGV);
my $markdown = ($ARGV[0] eq "-c")? 1 : 0;

if($markdown){
	shift(@ARGV);
	print "| FName | IntNum | TotLen | LenAvg | LenStdev | LenMedian | SmallestL | LargestL |\n";
	print "| :---- | -----: | -----: | -----: | -------: | --------: | --------: | -------: |\n";
}else{
	print "FName\tIntNum\tTotLen\tLenAvg\tLenStdev\tLenMedian\tSmallestL\tLargestL\n";
}
foreach my $f (@ARGV){
	open(IN, "< $f") || die "Could not open file $f!\n$usage";
	my @ints;
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		push(@ints, $segs[2] - $segs[1]);
	}
	close IN;
	my $int_num = scalar(@ints);
	my $int_tot = total_length(\@ints);
	my $int_avg = average(\@ints);
	my $int_stdev = standard_deviation(\@ints);
	my ($small, $large, $int_median) = median(\@ints);
	if($markdown){
		print "| $f | $int_num | $int_tot | $int_avg | $int_stdev | $int_median | $small | $large |\n";
	}else{
		print "$f\t$int_num\t$int_tot\t$int_avg\t$int_stdev\t$int_median\t$small\t$large\n";
	}
}

exit;

sub standard_deviation {
	my $array_ref = shift(@_);
#Prevent division by 0 error in case you get junk data
	my @numbers = @{$array_ref};
	return undef unless(scalar(@numbers));

# Step 1, find the mean of the numbers
	my $total1 = 0;
	foreach my $num (@numbers) {
	$total1 += $num;
	}
	my $mean1 = $total1 / (scalar @numbers);

# Step 2, find the mean of the squares of the differences
# between each number and the mean
	my $total2 = 0;
	foreach my $num (@numbers) {
	$total2 += ($mean1-$num)**2;
	}
	my $mean2 = $total2 / (scalar @numbers);

# Step 3, standard deviation is the square root of the
# above mean
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