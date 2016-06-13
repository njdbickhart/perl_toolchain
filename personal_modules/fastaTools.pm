#!/usr/bin/perl
# This is a collection of tools that I will use to encapsulate fasta file manipulation

package fastaTools;
use Mouse;
use namespace::autoclean;
use File::Basename;

has 'fasta' => (is => 'ro', isa => 'Str', required => 1);
has 'index' => (is => 'rw', isa => 'Str', builder => '_checkIdx');
has 'samInPath' => (is => 'rw', isa => 'Bool', builder => '_checkSamExe');

# A low-mem routine to estimate GC content from a given sequence
sub calcGCContent{
	my ($self, $chr, $start, $end) = @_;
		
	my $reference = $self->fasta;
	my $len = $end - $start;
	my $seq;
	if(!$self->samInPath){
		print STDERR "Error! Could not find samtools in path! Exiting!\n";
		exit;
	}
		
	open(my $IN, "samtools faidx $reference $chr:$start-$end |");
	my $h = <$IN>;
	print STDERR "$h";
	my $totGC = 0;
	while(my $line = <$IN>){
		chomp $line;
		my $gc = ($line =~ tr/GC/GC/);
		$totGC += $gc;
	}
	close $IN;
	return $totGC / $len;
}

# Takes single chr coordinate and returns just the sequence
sub retrieveSequence{
	my ($self, $chr, $start, $end) = @_;
	
	my $reference = $self->fasta;
	my $seq;
	if(!$self->samInPath){
		print STDERR "Error! Could not find samtools in path! Exiting!\n";
		exit;
	}
	
	open(my $IN, "samtools faidx $reference $chr:$start-$end |");
	my $h = <$IN>;
	while(my $line = <$IN>){
		chomp $line;
		$seq .= $line;
	}
	close $IN;
	return $seq;
}

sub _checkIdx{
	my ($self) = @_;
	my $supposedIdx = $self->fasta . ".fai";
	
	if(! -s $supposedIdx){
		if($self->samInPath){
			print STDERR "Could not find fasta index! Building now...\n";
			system("samtools faidx " . $self->fasta);
			return $self->fasta . ".fai";
		}else{
			print STDERR "Error! Could not find samtools in path! exiting!\n";
			exit;
		}
	}
	return $supposedIdx;		
}

sub _checkSamExe{
	my $found = 0;
	foreach my $path (split(/:/, $ENV{PATH})) {
		if( -f "$path/samtools") {
			$found = 1;
			last;
		}
	}
	return $found;
}

__PACKAGE__->meta->make_immutable;

1;