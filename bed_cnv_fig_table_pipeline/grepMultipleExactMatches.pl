#!/usr/bin/perl
# This script is designed to take a newline delimited list file or STDIN and to filter a file
# based on the input list

use strict;
use Getopt::Std;
use FileHandle;

my %opts;
my $usage = "perl $0 [options] -g <input list or STDIN> -c <column to search> -f <input file to search>
	This script is designed to select only exact matches from a tab file by searching distinct fields.
REQUIRED:
	-g	Input newline delimited entry file or \"stdin\" for STDIN
	-c	The column(s) to search from the [-f] file for exact matches (comma delimited numbers, zero-based)
	-f	The input file to search 
	
OPTIONAL:
	-o	Print data to this file [Default is to print to STDOUT]
	-x	Print only this/these columns from identified rows (comma delimited) [Default: print whole line]
	-e	Comment prefix. Exclude lines that begin with this character
\n";

getopt('gcfoxe', \%opts);
my $Obool = (defined($opts{'o'}))? 1 : 0;
my $Xbool = (defined($opts{'x'}))? 1 : 0;
my $Ebool = (defined($opts{'e'}))? 1 : 0;

unless(defined($opts{'g'}) && defined($opts{'c'}) && defined($opts{'f'})){
	print $usage;
	exit;
}

my @Scol = split(/,/, $opts{'c'});
my @Xcol; # columns to select from search file
if($Xbool){
	@Xcol = split(/,/, $opts{'x'});
}

# Determine input list type
my $il; # Input list filehandle
if($opts{'g'} eq "stdin"){
	$il = *STDIN;
}else{
	$il = FileHandle->new();
	$il->open("< $opts{g}") || die "Could not open input list: $opts{g}!\n$usage";
}

my %inputSet;
# Read in input list and store in set
while(my $line = <$il>){
	chomp $line;
	$inputSet{$line} = 1;
}

close $il;

# If we have output files specified, set these
my $out; # output filehandle
if($Obool){
	$out = FileHandle->new();
	$out->open("> $opts{o}");
}else{
	$out = *STDOUT;
}

# Now to loop through the input search file for the patterns we want
open(my $IN, "< $opts{f}") || die "Could not open input search file: $opts{f}!\n$usage";
while(my $line = <$IN>){
	if($Ebool){
		my $comment = $opts{'e'};
		if($line =~ /^$comment/){next;}
	}

	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line); 
	
	my $keep = 0; # boolean to print line
	foreach my $s (@Scol){
		if(exists($inputSet{$segs[$s]})){
			$keep = 1;
		}
	}
	
	# We found a hit! Now to process it
	if($keep){
		if($Xbool){
			my @retain;
			foreach my $x (@Xcol){
				push(@retain, $segs[$x]);
			}
			print {$out} join("\t", @retain) . "\n";
		}else{
			print {$out} "$line\n";
		}
	}
}

close $IN;
close $out;
exit;