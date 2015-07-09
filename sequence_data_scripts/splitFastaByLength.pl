#!/usr/bin/perl
# this is a one-shot script designed to create split fasta entries into sub 1 mb chunks

use strict;
use Getopt::Std;

my $usage = "perl $0 -f <input fasta> -o <output fasta file> \n";
my %opts;
getopt('fo', \%opts);

unless(defined($opts{'f'})){
	print $usage;
	exit;
}

unless( -s "$opts{f}.fai"){
	print STDERR "Error! Please index fasta!\n";
	exit;
}

# Get chr endpoints
open(my $IN, "< $opts{f}.fai") || die "Could not open fasta fai file!\n";
my %chrends;
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	
	$chrends{$segs[0]} = $segs[1];

}
close $IN;

open(my $OUT, "> $opts{o}");
open(my $IN, "< $opts{f}");
my $h = <$IN>;
chomp $h;
my $outh = $h;
my $ctr = 1;
if(NeedsCounter($h, \%chrends)){
	$outh = "$h.$ctr";
	$ctr++;
}
print $OUT "$outh\n";
my $bcounter = 0;
while(my $line = <$IN>){
	chomp($line);
	if($line =~ />/){
		$h = $line; 
		$outh = $h;
		$ctr = 1;
		if(NeedsCounter($h, \%chrends)){
			$outh = "$h.$ctr";
			$ctr++;
		}
		print $OUT "$outh\n";
		$bcounter = 0;
	}elsif($bcounter + length($line) >= 990000){
		print $OUT "$line\n";
		$outh = "$h.$ctr";
		$ctr++;
		print $OUT "$outh\n";
		$bcounter = 0;
	}else{
		print $OUT "$line\n";
		$bcounter += length($line);
	}
	
}
close $OUT;
close $IN;

exit;

sub NeedsCounter{
	my ($faheader, $key) = @_;
	my ($head) = $faheader =~ m/>(.+)/;
	if($key->{$head} > 990000){
		return 1;
	}else{
		return 0;
	}
}