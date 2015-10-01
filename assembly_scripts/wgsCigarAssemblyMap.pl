#!/usr/bin/perl
# This script aligns WGS reads from a fastq onto assembled contig fastas 
# The goal is to assess and concatenate CIGAR scores so that a bedgraph can be made

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -r <contig reference fasta> -f <wgs read fastq> -o <output bam file> -b <output bed file>\n";

getopt('rfob', \%opts);

unless(defined($opts{'r'}) && defined($opts{'f'}) && defined($opts{'o'}) && defined($opts{'b'})){
	print $usage;
	exit;
}

# Create the worker class
my $worker = CigarCondenser->new('reffasta' => $opts{'r'}, 'infastq' => $opts{'f'}, 'outbam' => $opts{'o'});

# perform alignments - generate output bam file
$worker->runAlignment();

# process the alignments into data strings for each base of the contig
$worker->processAlignments();

# produce output array
my $data = $worker->concatenateOutput();
open(my $OUT, "> $opts{b}");
foreach my $row (@{$data}){
	print {$OUT} join("\t", @{$row}) . "\n";
}
close $OUT;

exit;

BEGIN{
package CigarCondenser;
use Mouse;
use namespace::autoclean;

has ['reffasta', 'infastq', 'outbam'] => (is => 'ro', isa => 'Str', required => 1);
# scheme for data: {contig} -> {position} -> {CIGAR string} = read count
has 'datastore' => (is => 'rw', isa => 'HashRef[Any]', default => sub{{}});
has 'contiglens' => (is => 'rw', isa => 'HashRef[Any]', default => sub{{}});

sub runAlignment{
	my ($self) = @_;
	
	my $ref = $self->reffasta();
	my $fastq = $self->infastq();
	my $bam = $self->outbam();
	
	# Check if required indicies are present
	unless( -s "$ref.bwt"){
		print STDERR "BWA indexing ref fasta: $ref\n";
		system("bwa index $ref");
	}
	
	# Determine contig lens
	unless( -s "$ref.fai" ){
		print STDERR "Samtools faidx missing, preparing...\n";
		system("samtools faidx $ref");
	}
	
	my %contiglens;
	open(my $IN, "< $ref.fai");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$contiglens{$segs[0]} = $segs[1];
	}
	close $IN;
	
	$self->contiglens(\%contiglens);
	
	# output only mapping alignments
	my $sam = "temp" . int(rand(10000000)) . ".sam";
	system("bwa mem $ref $fastq | samtools view -bS -F 4 - | samtools sort -o $bam -T $bam.pre - ");
	#system("bwa aln $ref $fastq > $sam");
	#system("samtools view -bS -F 4 $sam | samtools sort -o $bam -T $bam.pre - ");
	
	print STDERR "Alignments finished. Bam file is: $bam\n";
	system("samtools index $bam");
}

sub processAlignments{
	my ($self) = @_;
	
	my $bam = $self->outbam();
	
	my %data;
	my $lastcontig; 
	open(my $IN, "samtools view $bam | ");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my $prevpos = $segs[3];
		foreach my $cigar ($segs[5] =~ /(\d+)([MDISH])/g){
			my $cstr = $2;
			my $ccount = $1;
			# Loop to assign each read coordinate position a cigar value counter
			for(my $x = $prevpos; $x < $prevpos + $ccount; $x++){
				if(exists($data{$segs[2]}->{$x})){
					if(exists($data{$segs[2]}->{$x}->{$cstr})){
						$data{$segs[2]}->{$x}->{$cstr} = $data{$segs[2]}->{$x}->{$cstr} + 1;
					}else{
						$data{$segs[2]}->{$x}->{$cstr} = 1;
					}
				}else{
					$data{$segs[2]}->{$x}->{$cstr} = 1;
				}
			}
			$prevpos += $ccount;
		}
	}
	close $IN;
	$self->datastore(\%data);	
}

sub concatenateOutput{
	my ($self) = @_;
	
	my $data = $self->datastore();
	my $contiglens = $self->contiglens;
	
	# output = []->[contig, start, end, cigarstring]
	my @output;
	
	foreach my $contig (sort{$a cmp $b} keys(%{$data})){
		my @rows;
		if(exists($contiglens->{$contig})){
			# we need to start the output array
			if(scalar(@output) < 1){
				if(!exists($data->{$contig}->{1})){
					push(@output, [$contig, 1, 1, "NONE"]);
				}else{
					my $cratio = $self->calculateCigarRatio($data->{$contig}->{1});
					push(@output, [$contig, 1, 1, $cratio]);
				}
			}
		
			for(my $pos = 1; $pos < $contiglens->{$contig}; $pos++){
			#foreach my $pos (sort{$a <=> $b} keys(%{$data->{$contig}}){
				if(!exists($data->{$contig}->{$pos})){
					# unmappable region of the contig
					if($output[-1]->[0] eq $contig && $output[-1]->[3] eq "NONE"){
						$output[-1]->[2] = $pos;
					}else{
						push(@output, [$contig, $pos, $pos, "NONE"]);
					}
				}else{
					# mapped region, process as normal
					my $cratio = $self->calculateCigarRatio($data->{$contig}->{$pos});
					if($output[-1]->[0] eq $contig && $output[-1]->[3] eq $cratio){
						$output[-1]->[2] = $pos;
					}else{
						push(@output, [$contig, $pos, $pos, $cratio]);
					}
				}
			}
		}
	}
				
	return \@output;
}

# I'm going to try cigar ratios (ie number of reads rounded up to 1 sig digits) first for concatenation
sub calculateCigarRatio{
	my ($self, $cigarHash) = @_;
	
	my $totreads = 0;
	foreach my $ckey (keys(%{$cigarHash})){
		$totreads += $cigarHash->{$ckey};
	}
	
	# It's sprintf concatenation time!
	# This enables me to reduce the sig digits for the ratio calculation
	# Ratios will be converted to reduced resolution percentages
	my $cstr = '';
	foreach my $ckey (sort {$a cmp $b} keys(%{$cigarHash})){
		my $ratio = sprintf("%.1f", $cigarHash->{$ckey} / $totreads);
		$ratio = int($ratio * 100);
		$cstr .= sprintf("%d%s", $ratio, $ckey);
	}
	
	return $cstr;
}
		

__PACKAGE__->meta->make_immutable;

}
