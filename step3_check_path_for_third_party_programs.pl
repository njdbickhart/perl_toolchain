#!/usr/bin/perl
# This script is designed to check the user's PATH to see if required third party programs are installed. 
# If any prerequisites are found to be missing, the user will just be informed at this time, as the weblinks are likely to change for
# each installed program.

use strict;
use Getopt::Std;

my @programlist = ("bwa", "samtools", "bcftools", "mrsfast");

checkReqs(@programlist);
exit;


sub checkReqs{
	my (@progs) = @_;
	
	my @notfound;
	foreach my $p (@progs){
		my $found = 0;
		foreach my $path (split(/:/, $ENV{PATH})) {
			if( -f "$path/$p") {
				$found = 1;
				last;
			}
		}
		if(!$found){
			push(@notfound, $p);
		}
	}
	if(scalar(@notfound) > 0){
		print "Error! Could not find the following programs on your path:\n";
		foreach my $p (@notfound){
			printf("%-24s.......NOT FOUND!\n", $p);
		}
		print "Please install the above third-party programs before using the scripts.\n";
		exit;
	}else{
		print "Checked for the presence of the following third party programs:\n";
		foreach my $p (@progs){
			printf("%-24s.......OK!\n", $p);
		}
		print "Your installation looks good!\n";
	}
}