#!/usr/bin/perl
# This script automatically parses Steve's automated directory structure in order to generate spreadsheets for
# use in my bwa runner program. I've improved the interface for the program and it should run almost instantly.

use strict;
use Class::Struct;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -d <base directory to scan> -o <output spreadsheet file> -k <key conversion for IDs>\n
	-d	This is the base directory to search for Steve's files (ie. the flowcell name)
	-o	Allows the user to designate the name for the output spreadsheet file (format: fq1(tab)fq2(tab)animalname(tab))
	-k	A key file to convert sample names (format: steveid(tab)newid). To convert \"lane\" into a proper id\n";
	
getopt('dok', \%opts);

# Check for help invocation
if($opts{'h'}){
	print $usage;
	exit;
}

# Options D and O are required
unless(defined($opts{'d'}) && defined($opts{'o'})){
	print "Missing mandatory arguments 'd' or 'o'!\n";
	print $usage;
	exit;
}

struct( read_sort => {
	'fq1' => '$',
	'fq2' => '$',
	'fnum' => '$',
	'fc' => '$',
	'segs' => '@',
});


my %data_struct; # {an_name}->[]->read_sort

# If the "k" option is invoked, read in the conversion table
my %key_conv;
if(defined($opts{'k'})){
	open(IN, "< $opts{k}") || die "Could not open key file: $opts{k}!\n$usage";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$key_conv{$segs[0]} = $segs[1];
	}
	close IN;
}

# Retrieve all of the fq files in the directories
my @files = getFqs($opts{'d'}, $usage);
# Check if we have any files to parse at all!
if(scalar(@files) == 0){
	print "Did not identify any fq files in the flowcell directory: $opts{d}!\n$usage";
	exit;
}

# Process the files sequentially and add them to the appropriate containers
foreach my $line (@files){
	chomp $line;
	# replace double directory lines with single lines
	$line =~ s/\/\//\//g;
	my @s = split(/\//, $line);
	my @r = split(/[_\.]/, $s[-1]);
	
	# Key conversion if the "k" option was invoked
	my $an;
	if(exists($key_conv{$r[0]})){
		print "Converting $r[0] to: " . $key_conv{$r[0]} . "\n";
		$an = $key_conv{$r[0]};
	}else {
		$an = $r[0];
	}
	
	my $found = 0;
	if($r[3] eq "R1"){
		# The file segment indicated that this was read1 of the pair
		if(exists $data_struct{$an}){
			foreach my $row (@{$data_struct{$an}}){
				if($row->fnum() eq $r[4] && $row->fc() eq $s[5]){
					# The mate was already in the data structure, so add it to the existing entry
					$row->fq1($line);
					$found = 1;
				}
			}
		}
		unless($found){
			# We didn't find the mate in our structure, so create a new one
			push(@{$data_struct{$an}}, read_sort->new(
				'fq1' => $line,
				'fnum' => $r[4],
				'segs' => \@r,
				'fc' => $s[5],
			));
		}
	}else{
		# The file segment indicated that this was read2 of the pair
		if(exists $data_struct{$an}){
			foreach my $row (@{$data_struct{$an}}){
				if($row->fnum() eq $r[4] && $row->fc() eq $s[5]){
					$row->fq2($line);
					$found = 1;
				}
			}
		}
		unless($found){
			push(@{$data_struct{$an}}, read_sort->new(
				'fq2' => $line,
				'fnum' => $r[4],
				'segs' => \@r,
				'fc' => $s[5],
			));

		}
	}
}
close IN;

#my %data_struct; # {an_name}->[]->read_sort
# $data_struct{$an}
open(OUT, "> $opts{o}");
foreach my $an (keys(%data_struct)){
	foreach my $row ( @{ $data_struct{$an} } ){
		print OUT $row->fq1() . "\t" . $row->fq2 . "\t$an\n";
	}
}
close OUT;

print "Done loading files into spreadsheet: $opts{o}\n";
exit;


sub getFqs{
	my ($dir, $usage) = @_;
	my @storage;
	# Open the base directory
	opendir(DIR, $dir) || die "Could not open directory: $dir!\n$usage";
	while(my $sampdir = readdir(DIR)){
		# Skip the by products of the Unix directory structure
		if($sampdir eq '..' || $sampdir eq '.'){next;}
		# Open up the sample directories for processing
		opendir(SAMP, "$dir/$sampdir") if(-d "$dir/$sampdir");
		while(my $file = readdir(SAMP)){
			# We're only interested in the files
			if(!(-f "$dir/$sampdir/$file")){next;}
			# Check for any file that ends with a fastq.gz, fq.gz, fastq or fq
			if($file =~ m/.+(fastq|fq)+(.gz)*/){
				push(@storage, "$dir/$sampdir/$file");
			}
		}
		# Parity checking for the directory
		if(scalar(@storage) == 0){
			print STDERR "Did not find fq files within directory: $dir/$sampdir!\n";
		}
		closedir(SAMP);
	}
	closedir(DIR);
	return @storage;
}