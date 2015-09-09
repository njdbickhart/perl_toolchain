#!/usr/bin/perl
# This script subsections fastq entries, aligns them to a reference genome and attempts to identify consensus assembly regions

use strict;

BEGIN{
package unitigFastq;
use Mouse;
use namespace::autoclean;

has ['sequence', 'origname', 'quality', 'tempfq'] => (is => 'rw', isa => 'Str');
has 'unitigMappings' => (is => 'rw', isa => 'ArrayRef[Any]');

sub loadFastq{
	my ($self, $name, $seq, $plus, $qual) = @_;
	
	my ($origname) = $name =~ m/\@(.+)$/;
	$self->origname($origname);
	$self->sequence($seq);
	$self->quality($qual);
	
	$self->tempfq("temp.fq");
	open(my $OUT, "> temp.fq");
	
	my $max = length($seq);
	for(my $x = 0; $x < $max + 1000; $x += 1000){
		my $tempend = $x + 1000;
		my $subseq = substr($seq, $x, 1000);
		my $subqual = substr($seq, $x, 1000);
		my $tempname = "$origname.$x.$tempend";
		
		print {$OUT} "\@$tempname\n$subseq\n+\n$subqual\n";
	}
	close $OUT;
	
}

sub alignAndProcess{
	my ($self, $reference) = @_;
	
	open(my $IN, "bwa mem $reference temp.fq |"); 
	while(my $line = <$IN>){
		

__PACKAGE__->meta->make_immutable;

package unitigMapping;
use Mouse;
use namespace::autoclean;

# Name contains the following attribute: fastqname.start.end
has ['name', 'chr'] => (is => 'rw', isa => 'Str');
has ['ustart', 'uend', 'start', 'end', 'quality', 'flag'] => (is => 'rw', isa => 'Num');

sub addSamLine{
	my ($self, $line) = @_;
	
	my @segs = split(/\t/, $line);
	my @namesegs = split(/\./, $segs[0]);
	if(scalar(@namesegs) < 3){
		print STDERR "Error with $segs[0] read!\n";
		return -1;
	}
	$self->name($namesegs[0]);
	$self->ustart($namesegs[1]);
	$self->uend($namesegs[2]);
	
	$self->chr($segs[2]);
	$self->start($segs[3]);
	$self->end($segs[3] + length($segs[9]));
	$self->quality($segs[4]);
	$self->flag($segs[1]);
}

__PACKAGE__->meta->make_immutable;

}