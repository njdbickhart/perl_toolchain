#!/usr/bin/perl
# This script subsections fastq entries, aligns them to a reference genome and attempts to identify consensus assembly regions
# Update - 9/14/2015:	Modified program to print out statistics on assembled segments at the end prior to printout.
# Update - 9/15/2015:	Added fasta support

use strict;
use Getopt::Std;

my $usage = "perl $0 -f <unitig fastq/fasta> -r <bwa indexed ref genome> -o <output tab file>\n";
my %opts;
getopt('fro', \%opts);

unless(defined($opts{'f'}) && defined($opts{'r'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $worker = unitigFastq->new();

open(my $IN, "< $opts{f}") || die "Could not open input: $opts{f}!\n";
# Determine file format
my $line = <$IN>;
seek($IN, 0, 0);

if($line =~ /^>/){
	# Fasta format
	my ($name, $seq);
	while(my $line = <$IN>){
		chomp $line;
		$line =~ s/\r//g;
		if($line =~ /^>/ && length($seq) > 1){
			$worker->loadFasta($name, $seq);
			$seq = '';
			$name = $line;
		}elsif($line =~ /^>/ && length($seq) < 1){
			$name = $line;
		}else{
			$seq .= $line;
		}
	}
	# Catch end of file fasta data
	if(length($name) > 1 && length($seq) > 1){
		$worker->loadFasta($name, $seq);
	}
	print STDERR "loaded fasta file!\n";
}elsif($line =~ /^@/){
	# Fastq format

	while(my $name = <$IN>){
		my $seq = <$IN>;
		my $plus = <$IN>;
		my $qual = <$IN>;

		chomp($name, $seq, $qual);
		$name =~ s/\r//g;
		$seq =~ s/\r//g;
		$qual =~ s/\r//g;
		$worker->loadFastq($name, $seq, $qual);
	}
	print STDERR "loaded fastq file!\n";
}else{
	# Indeterminate format!
	print STDERR "Error! Could not determine format of input text file!\n";
	exit;
}
close $IN;

$worker->alignAndProcess($opts{'r'});

open(my $OUT, "> $opts{o}");
my @outputlines = $worker->produceOutput();

foreach my $row (@outputlines){
	print {$OUT} join("\t", @{$row}) . "\n";
}
close $OUT;
print STDERR "Output in: $opts{o}\n";

exit;

BEGIN{
package unitigFastq;
use Mouse;
use namespace::autoclean;

has ['sequence', 'origname', 'quality', 'tempfq'] => (is => 'rw', isa => 'Str');
has 'unitigMappings' => (is => 'rw', isa => 'ArrayRef[Any]');

sub loadFastq{
	my ($self, $name, $seq, $qual) = @_;
	
	my ($origname) = $name =~ m/\@(.+)$/;
	$self->origname($origname);
	$self->sequence($seq);
	$self->quality($qual);
	
	$self->tempfq("temp.fq");
	open(my $OUT, "> temp.fq");
	
	my $max = length($seq);
	for(my $x = 0; $x < $max; $x += 1000){
		my $e = 1000;
		if($x + 1000 > length($seq)){
			$e = length($seq) - $x - 1;
		}
		my $tempend = $x + $e;
		my $subseq = substr($seq, $x, $e);
		my $subqual = substr($qual, $x, $e);
		my $tempname = "$origname.$x.$tempend";
		
		print {$OUT} "\@$tempname\n$subseq\n+\n$subqual\n";
	}
	close $OUT;
	
}

sub loadFasta{
	my ($self, $name, $seq) = @_;
	
	my ($origname) = $name =~ m/\>(.+)$/;
	$self->origname($origname);
	$self->sequence($seq);
	$self->quality('I' x length($seq));
	
	$self->tempfq("temp.fq");
	open(my $OUT, "> temp.fq");
		
	my $max = length($seq);
	for(my $x = 0; $x < $max; $x += 1000){
		my $e = 1000;
		if($x + 1000 > length($seq)){
			$e = length($seq) - $x - 1;
		}
		my $tempend = $x + $e;
		my $subseq = substr($seq, $x, $e);
		my $subqual = 'I' x $e;
		my $tempname = "$origname.$x.$tempend";

		print {$OUT} "\@$tempname\n$subseq\n+\n$subqual\n";
	}
	close $OUT;
}
	

sub alignAndProcess{
	my ($self, $reference) = @_;
	
	my @mappings;
	open(my $IN, "bwa mem $reference temp.fq |"); 
	while(my $line = <$IN>){
		if($line =~ /^\@/){next;}
		
		chomp $line;
		my $tempmap = unitigMapping->new();
		$tempmap->addSamLine($line);
		push(@mappings, $tempmap);
	}
	close $IN;
	
	@mappings = sort {$a->ustart <=> $b->ustart} @mappings;
	
	$self->unitigMappings(\@mappings);
}

sub produceOutput{
	my ($self) = @_;
	
	# name	ustart	uend	chr	start	end	qualitystring
	my @outlines;
	my $prevchr = "NA"; my $prevustart = 200000000; my $prevcstart = 200000000; my @qualities;
	my ($prevuend, $prevcend);
	foreach my $map (@{$self->unitigMappings}){
		if($map->chr() ne $prevchr && $prevchr ne "NA"){
			if($prevustart != 200000000 && $prevcstart != 200000000){
				push(@outlines, [$self->origname, $prevustart, $prevuend, $prevchr, $prevcstart, $prevcend, join(";", @qualities)]);
			}
			@qualities = ();
			push(@qualities, $map->quality);
			$prevchr = $map->chr();
			$prevustart = $map->ustart;
			$prevuend = $map->uend;
			$prevcstart = $map->start;
			$prevcend = $map->end;
		}elsif($prevchr eq "NA"){
			$prevchr = $map->chr();
			$prevustart = $map->ustart;
			$prevuend = $map->uend;
			$prevcstart = $map->start;
			$prevcend = $map->end;
			push(@qualities, $map->quality);
		}else{
			$prevuend = $map->uend;
			$prevcend = $map->end;
			push(@qualities, $map->quality);
		}
	}
	push(@outlines, [$self->origname, $prevustart, $prevuend, $prevchr, $prevcstart, $prevcend, join(";", @qualities)]);
	
	$self->printLongestChrAligns(\@outlines);
	
	return @outlines;
}

sub printLongestChrAligns{
	my ($self, $outlines) = @_;
	
	my %chrmins; my %chrmaxs; my %chrlens;
	foreach my $line (@{$outlines}){
		my $min = ($line->[4] < $line->[5])? $line->[4] : $line->[5];
		my $max = ($line->[4] > $line->[5])? $line->[4] : $line->[5];
		
		if(!exists($chrmins{$line->[3]})){
			$chrmins{$line->[3]} = $min;
			$chrmaxs{$line->[3]} = $max;
			$chrlens{$line->[3]} = $max - $min;
		}else{
			if($chrmins{$line->[3]} > $min){
				$chrmins{$line->[3]} = $min;
			}
			if($chrmaxs{$line->[3]} < $max){
				$chrmaxs{$line->[3]} = $max;
			}
			$chrlens{$line->[3]} = $chrmaxs{$line->[3]} - $chrmins{$line->[3]};
		}
	}	
	
	print "Longest aligments:\tchr\tstart\tend\tlength\n";
	foreach my $chrs (sort {$chrlens{$b} <=> $chrlens{$a}} keys(%chrlens)){
		print "\t\t\t$chrs\t" . $chrmins{$chrs} . "\t" . $chrmaxs{$chrs} . "\t" . $chrlens{$chrs} . "\n";
	}
}

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