#!/usr/bin/perl
# This script is designed to identify exons that span across contig alignments
# It should be used in conjunction with bed_cnv_fig_table_pipeline/generateExonFastaFromUCSCTable.pl

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -b <input bam> -c <[optional] conversion of IDs to genes> -o <output file>\n";
getopt('bco', \%opts);

unless(defined($opts{'b'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $convert = 0;
my %converter;
if(defined($opts{'c'})){
	$convert = 1;
	open(my $IN, "< $opts{c}") || die "Could not open conversion file!\n";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$segs[0] =~ s/\.\d{1,3}//;
		$converter{$segs[0]} = $segs[1];
	}
	close $IN;
}

# {gene} -> {contig} -> [] exon numbers
my %storage;
open(my $IN, "samtools view $opts{b} | ");
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	my($geneid, $exonnum) = $segs[0] =~ m/(.+)_(\d+)/;
	if($segs[2] =~ /\*/){
		$segs[2] = "unmapped";
	}
	if($convert){
		$geneid = $converter{$geneid};
	}
	push(@{$storage{$geneid}->{$segs[2]}}, $exonnum);
}
close $IN;

open(my $OUT, "> $opts{o}");
foreach my $genes (sort{scalar(@{$storage{$b}}) <=> scalar(@{$storage{$a}})} keys(%storage)){
	# Sort by number of different scaffold spans!
	print $OUT "$genes";
	foreach my $contigs (keys(%{$storage{$genes}})){
		my $exonstr = join(",", @{$storage{$genes}->{$contigs}});
		print $OUT "\t$contigs\t$exonstr";
	}
	print $OUT "\n";
}
close $OUT;
exit;