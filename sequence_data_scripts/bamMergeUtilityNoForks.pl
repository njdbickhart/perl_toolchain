#!/usr/bin/perl
# This script is designed to automate the merger of multiple bam files with read groups using samtools "Merge"
# It will also check for bam indexing and perform that first if it is missing
# This version does not use the forks::super module so that it can be compatible with the Ceres cluster

use strict;
use File::Basename;
use Getopt::Std;

my $usage = "Usage: perl $0 -i <comma delimited list of bam files to merge> [OR] -d <directory containing bam files to merge> -o <output merged bam name>
Select either:
	-i	A comma delimited list of bam file paths to merge
	OR
	-d	A directory containing bam files that will all be merged together
	
AND
	-o	The output bam file name (complete with path to the file)
	-r	FLAG: remove all merged bams if the merged bam is created and is not empty\n";

my %opts;
getopt('idon', \%opts);

if(defined($opts{'h'})){
	print $usage;
	exit;
}

unless((defined($opts{'i'}) || defined($opts{'d'})) && defined($opts{'o'})){
	print "Missing mandatory options!\n$usage";
	exit;
}

print "When life gives you John Coles, make Lemonade\n";

my $threads = 1;

# Fill required variables and flags
my @files;
my $remove = 0;

if(defined($opts{'i'})){
	chomp $opts{'i'};
	@files = split(/,/, $opts{'i'});
}elsif(defined($opts{'d'})){
	@files = `ls $opts{d}/*.bam`;
	chomp(@files);
}

if($opts{'r'}){
	$remove = 1;
}

# Check if bams need to be indexed
my @needsindexing;
foreach my $f (@files){
	unless(-s "$f.bai"){
		push(@needsindexing, $f);
	}
}

if(scalar(@needsindexing) > 0){
	print "Indexing bam files.\n";
	foreach my $f (@needsindexing){
		samIndexer($f);
	}
}
waitall();

# Create the header for the bam file
my $header = "header.temp.sam";
headerCreation(\@files, $header);

samMerge(\@files, $header, $opts{'o'}, $remove);

exit;

sub samMerge{
	my ($file_array, $header, $output, $remove) = @_;
	
	my $filestr = join(" ", @{$file_array});
	print "samtools merge -h $header $output $filestr\n";
	system("samtools merge -h $header $output $filestr");
	print "samtools index $output\n";
	system("samtools index $output");
	
	if($remove && -s $output){
		print "Cleaning up bam files...\n";
		foreach my $f (@{$file_array}){
			system("rm $f");
		}
	}
	
	if(-s $output){
		system("rm $header");
	}
}

sub headerCreation{
	my ($file_array, $header_name) = @_;
	my %order = ('@HD' => 0, '@SQ' => 1, '@RG' => 2, '@PG' => 3);
	my @store; # [0,1,2,3] -> []
	my $begin = 1;
	
	foreach my $f (@{$file_array}){
		unless(-s $f){
			print "Could not find file: $f for samtools header creation!\n";
		}
		
		open(IN, "samtools view -H $f | ");
		while(my $line = <IN>){
			chomp $line;
			my @segs = split(/\t/, $line);
			
			if(!exists($order{$segs[0]})){next;}
			if($begin && $segs[0] ne '@RG'){
				my $num = $order{$segs[0]};
				push(@{$store[$num]}, $line);
			}elsif($segs[0] eq '@RG'){
				push(@{$store[2]}, $line);
			}
		}
		close IN;
		$begin = 0;
	}
	
	open(OUT, "> $header_name");
	for(my $x = 0; $x < scalar(@store); $x++){
		foreach my $row (@{$store[$x]}){
			print OUT "$row\n";
		}
	}
	close OUT;
}
sub samIndexer{
	my ($file) = @_;
	print "Could not find bam index for: $file. Indexing now...\n";
	system("samtools index $file");
}