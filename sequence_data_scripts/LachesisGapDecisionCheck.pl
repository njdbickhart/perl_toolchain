#!/usr/bin/perl
# This script takes regions of known contig/scaffold splits and intersects them with 
# perpendicular datasets to make a decision for orientation

use strict;
use GetOpt::Std;

my %opts;
getopt('cgbfuo', \%opts);
my $usage = "perl $0 -c <contig list coords> -g <gap coords> -b <BND interchrom bedpe> -f <filled gaps> -u <unfilled gaps> -o <output base>\n";

unless(defined($opts{'c'}) && defined($opts{'b'}) && defined($opts{'f'}) && defined($opts{'u'}) && defined($opts{'o'})){
	 print $usage;
	 exit;
}

my $worker = eventWorkhorse->new('filename' => $opts{'g'});
$worker->loadGapInfo();
print STDERR "Loaded gap site file\n";

$worker->associateWithSites($opts{'b'}, $opts{'f'}, $opts{'u'});
print STDERR "Finished association\n";

$worker->printOutGapSites($opts{'o'});
print STDERR "Printed output to file: $opts{o}.gap.tab\n";
exit;

BEGIN{
package eventWorkhorse;
use Mouse;
use namespace::autoclean;

our @FILETYPES = ('FILLED', 'UNFILLED', 'BND');

# hash ref is ordered by name 
has 'data' => (is => 'rw', isa => 'HashRef[event]');
has 'filename' => (is => 'ro', isa => 'Str');
# key-> file type -> name -> events array
has 'associations' => (is => 'rw', isa => 'HashRef[HashRef[event]]');

sub loadGapInfo{
	my ($self) = @_;
	
	my %data;
	open(my $IN, "< $self->filename") || die "Could not open input file!\n";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$data{$segs[3]} = event->new('chr'=> $segs[0], 
			'start' => $segs[1],
			'end' => $segs[2],
			'name' => $segs[3]);
	}
	close $IN;
}

sub printOutGapSites{
	my ($self, $output) = @_;
	
	my $data = $self->data;
	my $assoc = $self->associations;
	open(my $OUT, "> $output.gap.tab");
	foreach my $name (keys(%{$data})){
		# Check what category this event fills in with
		my @types;
		if(exists($assoc->{'FILLED'}->{$name})){
			push(@types, 'FILLED');
		}elsif(exists($assoc->{'UNFILLED'}->{$name})){
			push(@types, 'UNFILLED');
		}elsif(exists($assoc->{'BND'}->{$name})){
			push(@types, $assoc->{'BND'}->{$name}->name);
		}
		
		my $ref = $data->{$name}->getOutArray;
		push(@{$ref}, join("/", @types));
		print {$OUT} join("\t", @{$ref}) . "\n";
	}
	close $OUT;
}

sub associateWithSites{
	my ($self, $bnd, $filled, $unfilled) = @_;
	my $file = $self->filename;
	
	# Easy filled sites
	my %assoc;
	open( my $IN, "intersectBed -a $file -b $filled |");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$assoc{"FILLED"}->{$segs[3]} = event->new('chr'=> $segs[0], 
			'start' => $segs[1],
			'end' => $segs[2],
			'name' => $segs[3]);
	}
	close $IN;
	
	# Easy unfilled sites
	open($IN, "intersectBed -a $file -b $unfilled |");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$assoc{"UNFILLED"}->{$segs[3]} = event->new('chr'=> $segs[0], 
			'start' => $segs[1],
			'end' => $segs[2],
			'name' => $segs[3]);
	}
	close $IN;
	
	# complex bnd events
	open($IN, "intersectBed -a $file -b $bnd -wb |");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my $ename;
		my $fname = $segs[4]; my $sname = $segs[7];
		
		# scaffold/contig name + orientation
		if($fname eq $sname){
			$ename = "$fname-$segs[12];$segs[13]";
		}else{
			$ename = "$fname;$sname-$segs[12];$segs[13]";
		}
		$assoc{"BND"}->{$segs[3]} = event->new('chr'=> $segs[0], 
			'start' => $segs[1],
			'end' => $segs[2],
			'name' => $ename);
	}
	close $IN;
	$self->associations(\%assoc);
}
	

__PACKAGE__->meta->make_immutable;

package event;
use Mouse;
use namespace::autoclean;

has ['chr','name'] => (is => 'ro', isa => 'Str');
has ['start', 'end'] => (is => 'ro', isa => 'Int');

sub getOutArray{
	my ($self) = @_;
	return [$self->chr, $self->start, $self->end, $self->name];
}

__PACKAGE__->meta->make_immutable;
}