#!/usr/bin/perl
# this is a one-shot script designed to create split fasta entries from a bed file with breakpoint coordinates

use strict;
use Getopt::Std;

my $usage = "perl $0 -f <indexed fasta file> -o <output fasta file> -b <bed file with breakpoints>\n";
my %opts;
getopt('fob', \%opts);

unless(defined($opts{'f'})){
	print $usage;
	exit;
}

unless( -s "$opts{f}.fai"){
	print STDERR "Error! Please index fasta!\n";
	exit;
}

my %coords; # {ctg} -> [] ->[start, end]
my %chrends; # {ctg} = endcoord
# store the bed coordinates
open(my $IN, "< $opts{b}") || die "Could not open bed file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	push(@{$coords{$segs[0]}}, [$segs[1], $segs[2]]);
}
close $IN;

# Get chr endpoints
open(my $IN, "< $opts{f}.fai") || die "Could not open fasta fai file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	if(exists($coords{$segs[0]})){
		$chrends{$segs[0]} = $segs[1];
	}
}
close $IN;

open(my $OUT, "> $opts{o}");
# Work on split contigs first
foreach  my $ctg (keys(%coords)){
	my @breaks = @{$coords{$ctg}};
	
	if(scalar(@breaks) == 1){
		my $start = $breaks[0]->[0];
		my $end = $breaks[0]->[1];
		open(my $IN, "samtools faidx $opts{f} $ctg\:1-$start |");
		my $h = <$IN>;
		print $OUT ">$ctg.1\n";
		while(my $line = <$IN>){
			print $OUT $line;
		}
		close $IN;
		
		my $chrend = $chrends{$ctg};
		open(my $IN, "samtools faidx $opts{f} $ctg\:$end-$chrend |");
		my $h = <$IN>;
		print $OUT ">$ctg.2\n";
		while(my $line = <$IN>){
			print $OUT $line;
		}
		close $IN;
	}else{
		my $prev = 1;
		my $ctr = 0;
		for(my $x = 0; $x <= scalar(@breaks); $x++){
			my $end = $breaks[$x]->[0];
			$ctr++;
			open(my $IN, "samtools faidx $opts{f} $ctg\:$prev-$end |");
			my $h = <$IN>;
			print $OUT ">$ctg.$ctr\n";
			while(my $line = <$IN>){
				print $OUT $line;
			}
			close $IN;
			$prev = $breaks[$x]->[1];
		}
		
		my $chrend = $chrends{$ctg};
		$ctr++;
		open(my $IN, "samtools faidx $opts{f} $ctg\:$prev-$chrend |");
		my $h = <$IN>;
		print $OUT ">$ctg.$ctr\n";
		while(my $line = <$IN>){
			print $OUT $line;
		}
		close $IN;
	}
}

print STDERR "Working on non-altered fasta headers\n";
# Now, open up the whole fasta file and only print it out if the contig didn't need to be split
open(my $IN, "< $opts{f}");
my $print = 0;
while(my $line = <$IN>){
	chomp $line;
	if($line =~ /^>/){
		my ($ctg) = $line =~ />(.+)/;
		if(!exists($coords{$ctg})){
			$print = 1;
			print $OUT "$line\n";
		}else{
			$print = 0;
		}
	}elsif($print){
		print $OUT "$line\n";
	}
}
close $IN;
close $OUT;

exit;