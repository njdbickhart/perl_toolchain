#!/usr/bin/perl
# This script sorts a bed file by chromosome, then by the start coordinate
# Takes stdin input

use strict;

my %chrsholder;
while(my $line = <STDIN>){
	chomp $line;
	my @s = split(/\t/, $line);
	push(@{$chrsholder{$s[0]}}, [@s]);
}

foreach my $chr (sort {
			my ($c1, $x) = $a =~ /(chr)*(.+)/;
			my ($c2, $y) = $b =~ /(chr)*(.+)/;
			if($x eq "X"){
				$x = 500;
			}elsif($x eq "Y"){
				$x = 501;
			}elsif($x eq "M" || $x eq "MT"){
				$x = 502;
			}
			
			if($y eq "X"){
				$y = 500;
			}elsif($y eq "Y"){
				$y = 501;
			}elsif($y eq "M" || $x eq "MT"){
				$y = 502;
			}
			$x <=> $y}
			keys(%chrsholder)){
	foreach my $row (sort {$a->[1] <=> $b->[1]} @{$chrsholder{$chr}}){
		print join("\t", @{$row}) . "\n";
	}	
}

exit;
