{
package genomic_coord;
use strict;
use Mouse;
use namespace::autoclean;

has 'chr' => ( is => 'rw', isa => 'Str',);
has 'start' => ( is => 'rw', isa => 'Int',);
has 'end' => ( is => 'rw', isa => 'Int',);
has 'name' => ( is => 'rw', isa => 'Str',);

__PACKAGE__->meta->make_immutable;
}
{
package cnvr_coord;
use strict;
use Mouse;
use namespace::autoclean;
extends 'genomic_coord';

#has 'anname' => (is => 'rw', isa => 'ArrayRef[Str]', builder => '_split_aname', lazy => 1);
#has 'num_unique_ans' => (is => 'rw', isa => 'Int', builder => '_count_unique', lazy => 1);
#has 'type' => (is => 'rw', isa => 'ArrayRef[Str]', builder => '_split_tname', lazy => 1);
has 'loci_hit' => (is => 'rw', isa => 'HashRef[ArrayRef[overlap_values]]');

sub split_aname{
	my $self = shift(@_);
	my @ret;
	my @s = split(/;/, $self->name());
	foreach my $r (@ret){
		my @t = split(/_/, $r);
		push(@ret, $t[0]);
	}
	return \@ret;
}
sub split_tname{
	my $self = shift(@_);
	my @ret;
	my @s = split(/;/, $self->name());
	foreach my $r (@ret){
		my @t = split(/_/, $r);
		push(@ret, $t[1]);
	}
	return \@ret;
}
sub count_unique{
	my $self = shift(@_);
	my %h;
	foreach my $v (@{$self->anname()}){
		$h{$v} = 1;
	}
	return scalar(keys(%h));
}

__PACKAGE__->meta->make_immutable;
}
{
package annotation_file;
use strict;
use Mouse;
use File::Basename;
use lib dirname(__FILE__);
use kentBinTools;
use namespace::autoclean;

has 'file' => ( is => 'rw', isa => 'Str', predicate => 'has_file',);
has 'bed' => (is => 'rw', isa => 'HashRef[HashRef[ArrayRef[genomic_coord]]]', predicate => 'has_bed',); # {chr}->{bin}->[0]->genomic_coord
has 'lines' => (is => 'rw', isa => 'Int',);

sub load_file {
	my $self = shift;
	if($self->has_file()){
		open(IN, "< $self->file()") || die "Could not open $self->file()!\n";
		while (my $line = <IN>){
			chomp  $line;
			$line =~ s/\r//g;
			if($line eq '' || $line =~ /^\s+$/){next;}
			my @segs = split(/\t/, $line);
			if(!(scalar(@segs) >= 4)){
				print "Error! Line has less than four categories! $line\n";
			}
			
			my $entry = genomic_coord->new('chr' => $segs[0], 'start' => $segs[1], 'end' => $segs[2], 'name' => $segs[3]);
			my $bin = kent_bin_tools::getbin($segs[1], $segs[2]);
			if(!($self->has_bed())){
				my %h;
				push(@{$h{$segs[0]}->{$bin}}, $entry);
				$self->bed(\%h);
			}else{
				my %h = $self->bed();
				push(@{$h{$segs[0]}->{$bin}}, $entry);
				$self->bed(\%h);
			}
			$self->lines($self->lines() + 1);
		}
		close IN;
	}else{
		print "Error! Did not set file attribute in annotation_file class! $self->file()\n";
		exit;
	}
	#print "Loaded file: $self->file()\n";
}
# Bool: does the input information intersect with the stored data
sub data_intersects {
	my ($self, $chr, $start, $end) = @_;
	unless($self->has_bed()){
		print "Error! Called data_intersects before loading data! $self->file()\n";
		exit;
	}
	if(!(exists($self->bed()->{$chr}))){
		return 0;
	}
	
	my @bins = kent_bin_tools::searchbins($start, $end);
	foreach my $b (@bins){
		if(!(exists($self->bed()->{$chr}->{$b}))){
			next;
		}
		foreach my $gc (@{$self->bed()->{$chr}->{$b}}){
			if(kent_bin_tools::overlap($gc->start(), $gc->end(), $start, $end) > 1){
				return 1;
			}
		}
	}
	return 0;
}

# Str, float, float: returns the name of the entry that overlapped, as well as the overlapping percentages
sub retrieve_info_intersection {
	my($self, $chr, $start, $end) = @_;
	my($ovlp_name, $a_perc, $b_perc);
	unless($self->has_bed()){
		print "Error! Called data_intersects before loading data! $self->file()\n";
		exit;
	}
	if(!(exists($self->bed()->{$chr}))){
		return "NA", 0, 0;
	}
	my @bins = kent_bin_tools::searchbins($start, $end);
	foreach my $b (@bins){
		if(!(exists($self->bed()->{$chr}->{$b}))){
			next;
		}
		foreach my $gc (@{$self->bed()->{$chr}->{$b}}){
			my $ov = kent_bin_tools::overlap($gc->start(), $gc->end(), $start, $end);
			if($ov > 1){
				my $alen = $end - $start;
				my $blen = $gc->end() - $gc->start();
				$a_perc = $ov / $alen;
				$b_perc = $ov / $blen;
				$ovlp_name = $gc->name();
				return $ovlp_name, $a_perc, $b_perc;
			}
		}
	}
	return "NA", 0, 0;
}
sub DEMOLISH {
	#print "Clearing file from memory space...\n";
}
__PACKAGE__->meta->make_immutable;
}
{
package overlap_values;
use Mouse;
use namespace::autoclean;

has 'ovlp_name' => (is => 'rw', isa => 'Str');
has 'cnvr_perc' => (is => 'rw', isa => 'Num');
has 'loci_perc' => (is => 'rw', isa => 'Num');

__PACKAGE__->meta->make_immutable;
}
{
package bin_bed_container;
use Mouse;
use File::Basename;
use lib dirname(__FILE__);
use kentBinTools;
use namespace::autoclean;

has 'bed' => (is => 'rw', isa => 'HashRef[HashRef[ArrayRef[Any]]]', default => sub{ {} },); # {chr}->{bin}->[0]->genomic_coord

sub add_bed{
	my ($self, $chr, $start, $end, $name) = @_;
	my $gcoord = genomic_coord->new(
		'chr' => $chr,
		'start' => $start,
		'end' => $end,
		'name' => $name,);
	my $bin = kent_bin_tools::getbin($start, $end);
	push(@{$self->bed->{$chr}->{$bin}}, $gcoord);
}



__PACKAGE__->meta->make_immutable;
}
{
package unbin_bed_container;
use Mouse;
use File::Basename;
use namespace::autoclean;

has 'bed' => (is => 'rw', isa => 'HashRef[ArrayRef[Any]]]', default => sub{ {} }, ); # {chr}->[0]->genomic_coord

sub add_bed{
	my ($self, $chr, $start, $end, $name) = @_;
	my $gcoord = genomic_coord->new(
		'chr' => $chr,
		'start' => $start,
		'end' => $end,
		'name' => $name,);
	push(@{$self->bed->{$chr}}, $gcoord);
}

sub get_ascend_sortedbed{
	my $self = shift;
	my $working = $self->bed();
	
	my @sorted;
	foreach my $chr (sort {
		my ($x) = $a =~ /chr(.+)/;
		my ($y) = $b =~ /chr(.+)/;
		if($x eq "X"){
			$x = 500;
		}elsif($x eq "Y"){
			$x = 501;
		}
		if($y eq "X"){
			$y = 500;
		}elsif($y eq "Y"){
			$y = 501;
		}
		$x <=> $y}
		keys(%{$working})){
		
		foreach my $row (sort { $a->start() <=> $b->start()} @{$working->{$chr}}){
			push(@sorted, $row);
		}		
	}
	return @sorted;
}

sub get_descend_sortedbed{
	my $self = shift;
	my $working = $self->bed();
	
	my @sorted;
	foreach my $chr (sort {
		my ($x) = $a =~ /chr(.+)/;
		my ($y) = $b =~ /chr(.+)/;
		if($x eq "X"){
			$x = 500;
		}elsif($x eq "Y"){
			$x = 501;
		}
		if($y eq "X"){
			$y = 500;
		}elsif($y eq "Y"){
			$y = 501;
		}
		$y <=> $x}
		keys(%{$working})){
		
		foreach my $row (sort { $b->start() <=> $a->start()} @{$working->{$chr}}){
			push(@sorted, $row);
		}		
	}
	return @sorted;
}

__PACKAGE__->meta->make_immutable;
}
1;