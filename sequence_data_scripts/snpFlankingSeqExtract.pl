#!/usr/bin/perl

use strict;
use Getopt::Std;
use Class::Struct;

# Set "autoflush" option so that the data is automatically printed
$| = 1;

# A way to organize our variables
struct(bedcoord => {
'start' => '$',
'len' => '$',
});

my $usage = "$0 -i <input reference fasta> -b <input snp coordinate file> -w <bases upstream and downstream> -o <output fasta>\n";

my %opts;
getopt('ibwoh', \%opts);

if(defined($opts{h})){
	print $usage;
	exit;
}

# Exit unless we find all of the variables we need!
unless(defined($opts{i}) && defined($opts{b}) && defined($opts{w}) && defined($opts{o})){
	print $usage;
	exit;
}

my %bedcoords; # {chr} ->[] = bedcoord

# First, read in the SNP coordinates so that we can zip through the fasta file later
open(IN, "< $opts{b}") || die "Could not open snp coordinate file!\n $usage";
while(my $line = <IN>){
	chomp $line;
	$line =~ s/\r//g;
	my @segs = split(/\t/, $line);
	# The coordinates are stored in this structure for easy access later
	my $coord = bedcoord->new(
		'start' => ($segs[1] - $opts{w} > 0)? $segs[1] - $opts{w} : 0,
		'len' => ($opts{w} + $opts{w}),
		);
	push(@{$bedcoords{$segs[0]}}, $coord);
}
close IN;

print "Snp coordinates loaded... beginning reference fasta processing\n";

# Now we run through the fasta file and try to pull the flanking sequence from each SNP site
open(FA, "< $opts{i}") || die "Could not open reference fasta!\n $usage";
open(OUT, "> $opts{o}") || die "Could not open output file!\n $usage";
my ($chr, $seq);
while(my $line = <FA>){
	chomp $line;
	# Fasta format is not clearly delimited. We need to store the fasta header and 
	# process lines until we reach the next header!
	if($line =~ />/){
		if(!defined($chr)){
			# First line! Store the chr header
			($chr) = $line =~ />(.+)/;
			print "Working on\t$chr\r";
			next;
		}else{
			# We reached a new fasta header. Process the data we've stored in the meantime
			my @falines = associate_coord_fasta($chr, $seq, \%bedcoords);
			if(scalar(@falines) > 0){
				foreach my $f (@falines){
					print OUT $f;
				}
			}
			($chr) = $line =~ />(.+)/;
			print "Working on\t$chr\r";
			$seq = '';
		}
	}else{
		# When we're not in the fasta header, store the sequence for later retrieval
		$seq .= $line;
	}
}
if(length($seq) > 0){
	# The end of the file doesn't have a fasta header, so we must print out the sequence afterwards!
	my @falines = associate_coord_fasta($chr, $seq, \%bedcoords);
	foreach my $f (@falines){
		print OUT $f;
	}
}

print "\nOutput fasta sequence in $opts{o}\n";

exit;

# This subroutine processes a sequence string to generate the flanking sequence we need
sub associate_coord_fasta{
	my ($chr, $seq, $bedref) = @_;
	my @ret;
	foreach my $coord (@{$bedref->{$chr}}){
		my $end = $coord->start() + $coord->len();
		my $faname = ">" . $chr . ":" . $coord->start() . "-" . $end;
		my $subseq = substr($seq, $coord->start() + 1, $coord->len());
		push(@ret, "$faname\n$subseq\n");
	}
	return @ret;
}