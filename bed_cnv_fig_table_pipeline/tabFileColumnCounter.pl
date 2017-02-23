#!/usr/bin/perl
# This script is designed to process a tab delimited file to generate some basic statistics on column field count
# It only works on columns that have a limited range of possible entries (like enumerated types)
# 6/12/2015: added a STDIN functionality
# 4/18/2016: added feature to skip empty lines
# 12/12/2016: made it so that all single space characters are delimited by default
# 2/23/2017: added numerical counting feature -- not tested on table generation!!!

use strict;
use Getopt::Std;
use columnCounter;

my $usage = "perl $0 [options] -f <input file> -c <column number>
	This script is designed to process a tab delimited file to generate some basic statistics on column field count.
REQUIRED:
	-f	Input space/tab delimited files separated by commas OR \"stdin\" for STDIN. Beginning spaces are removed from each line. Trailing spaces are condensed
	-c	The column number (Zero based!) to count features
	
OPTIONAL:
	-o	Output file [Default option: STDOUT]
	-e	Comment prefix. Exclude lines that begin with this character [Default option: exclude nothing]
	-m	Markdown flag. Formats output into table format [Default option: tab delimited]
	-x	Numerical counter statistics (different logic, calculates stdev, avg, median, etc) [Default option: type counter]
\n";

my %opts;
getopt('fcoex', \%opts);

unless(defined($opts{'f'}) && defined($opts{'c'})){
	print $usage;
	exit;
}

my $mbool = (defined($opts{'m'}))? 1 : 0;

# Determine how many files we're working with here.
my @files = split(/,/, $opts{'f'});
if(scalar(@files) == 1){
	# Just a single comparison
	my $worker = columnCounter->new('colnum' => $opts{'c'}, 'mkdwn' => $mbool);
	if(defined($opts{'x'})){
		$worker->numeric(1);
	}
	
	if(defined($opts{'e'})){
		$worker->ignore($opts{'e'});
	}
	
	if(defined($opts{'o'})){
		$worker->output($opts{'o'});
	}
	$worker->readFile($opts{'f'});
	$worker->writeOut();
}else{
	# Now we have more to work with
	my $manager = SummaryManager->new('mkdwn' => $mbool);
	
	if(defined($opts{'o'})){
		$manager->output($opts{'o'});
	}
	for (my $x = 0; $x < scalar(@files); $x++){
		my $filename = "File" . ($x + 1);
		
		my $worker = columnCounter->new('colnum' => $opts{'c'}, 'mkdwn' => $mbool);
		
		if(defined($opts{'x'})){
			$worker->numeric(1);
		}
		
		if(defined($opts{'e'})){
			$worker->ignore($opts{'e'});
		}
		
		$worker->readFile($files[$x]);
		$manager->setTable($filename, $worker->createSumTable($filename));
	}
	
	if(defined($opts{'o'})){
		open(OUT, "> $opts{o}");
		for(my $x = 0; $x < scalar(@files); $x++){
			my $filename = "File" . ($x + 1);
			print OUT "$filename\:  $files[$x]\n";
		}
		print OUT "\n";
		close OUT;
	}else{
		for(my $x = 0; $x < scalar(@files); $x++){
			my $filename = "File" . ($x + 1);
			print "$filename\:  $files[$x]\n";
		}
		print "\n";
	}
	
	$manager->PrintResults();
}

exit;

