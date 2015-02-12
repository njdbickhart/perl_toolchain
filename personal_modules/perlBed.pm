#!/usr/bin/perl
# A collection of classes used to manipulate bed files in Perl

package BedCoord;
use strict;
use Mouse;
use namespace::autoclean;

has 'chr' => ( is => 'ro', isa => 'Str',);
has 'start' => ( is => 'ro', isa => 'Int',);
has 'end' => ( is => 'ro', isa => 'Int',);
# Package for when the order of values is important
has 'other' => ( is => 'rw', isa => 'ArrayRef[Any]', predicate => 'has_other',);
# Package for when the values need easy lookup
has 'values' => ( traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', predicate => 'has_values', 
	handles => {'get_value' => 'get', 'num_values' => 'count', 'key_value' => 'keys'});

# can accept a list of ordered keys if a hash value was generated for the bed
sub outString{
	my ($self, $keyorder) = @_;
	my $c = $self->chr;
	my $s = $self->start;
	my $e = $self->end;
	
	if($self->has_other){
		my $o = $self->other;
		return "$c\t$s\t$e\t" . join("\t", @{$o}) . "\n";
	}elsif($self->has_values){
		my @data;
		if(defined($keyorder)){
			# Sort by input key order
			foreach my $k (@{$keyorder}){
				push(@data, $self->get_value($k));
			}
		}else{
			# Just sorting by lexicographical order
			foreach my $k (sort {$a cmp $b} $self->keys){
				push(@data, $self->get_value($k));
			}
		}
		return "$c\t$s\t$e\t" . join("\t", @data) . "\n";
	}else{
		return "$c\t$s\t$e\n";
	}
}

__PACKAGE__->meta->make_immutable;

package BedContainer;
use strict;
use Mouse;
use kentBinTools;
use namespace::autoclean;

has 'bed' => (is => 'rw', isa => 'HashRef[Any]', predicate => 'has_bed',); # {chr}->{bin}->[0]->genomic_coord
has 'lines' => (is => 'rw', isa => 'Int', default => 0,);


# takes longer, but is sorted based on chromosome and then start coordinate position
sub getOrderedList{
	my ($self) = @_;
	my %chrsholder = %{$self->bed};
	my @return;
	foreach my $chr (sort {
			my ($x) = $a =~ /chr(.+)/;
			my ($y) = $b =~ /chr(.+)/;
			if($x eq "X"){
				$x = 500;
			}elsif($x eq "Y"){
				$x = 501;
			}elsif($x eq "M" || $x eq "MT"){
				$x = 502;
			}
			
			if($y eq "X"){
				$y = 500;
			}elsif($y eq "Y"){
				$y = 501;
			}elsif($y eq "M" || $x eq "MT"){
				$y = 502;
			}
			$x <=> $y}
			keys(%chrsholder)){
		my @beds;
		foreach my $b (keys %{$chrsholder{$chr}}){
			push(@beds, @{$chrsholder{$chr}->{$b}});
		}
		push(@return, (sort{$a->start <=> $b->start} @beds));
	}
	return @return;
}

# Just a straight data dump with no sorting
sub getUnorderedList{
	my ($self) = @_;
	my @beds;
	my $ref = $self->bed;
	
	foreach my $chr (keys(%{$ref})){
		foreach my $bin (keys(%{$ref->{$chr}})){
			foreach my $bed (@{$ref->{$chr}->{$bin}}){
				push(@beds, $bed);
			}
		}
	}
	return @beds;
}

# Used for merging containers together
sub duplicateContainer{
	my ($self, $container) = @_;
	
	unless($container->isa('BedContainer')){
		print STDERR "[perlBed] Error! Passed object is not a BedContainer!\n";
		return;
	}
	
	if($self->has_bed){
		# Duplicate this the hard way...
		my $ref = $self->bed;
		my $temp = $container->bed;
		
		foreach my $chr (keys %{$temp}){
			foreach my $bin (keys %{$temp->{$chr}}){
				foreach my $bed (@{$temp->{$chr}->{$bin}}){
					push(@{$ref->{$chr}->{$bin}}, $bed);
				}
			}
		}
		$self->bed($ref);
	}else{
		# Easy duplication
		$self->bed($container->bed);
	}
}

# Less efficient way of adding bed reference to container
# Only use if infrequent pushes to the container
sub addBed{
	my ($self, $bed) = @_;
	
	my $binner = kentBinTools->new();
	my $ref = $self->bed;
	
	unless($bed->isa('BedCoord')){
		print STDERR "[perlBed] Error! Attempted to add a non-BedCoord object to the container!\n";
		exit;
	}
	my $bin = $binner->getbin($bed->start, $bed->end);
	push(@{$ref->{$bed->chr}->{$bin}}, $bed);
	
	$self->bed($ref);
}
# More efficient way of adding beds to the container
sub addBedArray{
	my ($self, $bedarray) = @_;
	
	my $binner = kentBinTools->new();
	my $ref = $self->bed;
	foreach my $bed (@{$bedarray}){
		unless($bed->isa('BedCoord')){
			print STDERR "[perlBed] Error! Attempted to add a non-BedCoord object to the container!\n";
			exit;
		}
		my $bin = $binner->getbin($bed->start, $bed->end);
		push(@{$ref->{$bed->chr}->{$bin}}, $bed);
	}
	
	$self->bed($ref);
}

# Requires the file name
# Optional: an arrayref of column key values for hashing
sub loadFile {
	my ($self, $file, $colkeys) = @_;
	
	my $binner = kentBinTools->new();
	open(IN, "< $file") || die "[perlBed] Could not open $file!\n";
	while (my $line = <IN>){
		chomp  $line;
		$line =~ s/\r//g;
		if($line eq '' || $line =~ /^\s+$/){next;}
		my @segs = split(/\t/, $line);
		if(!(scalar(@segs) >= 3)){
			print STDERR "[perlBed] Error! Line has less than four categories! $line\n";
			exit;
		}
		
		my $segNum = scalar(@segs);
		my $entry;
		if(defined($colkeys) && scalar(@segs > 3)){
			my %h;
			for (my $x = 0; $x < scalar(@{$colkeys}); $x++){
				my $y = $x + 3;
				$h{$colkeys->[$x]} = $segs[$y];
			}
			$entry = BedCoord->new('chr' => $segs[0], 'start' => $segs[1], 'end' => $segs[2], 'values' => \%h);
		}elsif(!defined($colkeys) && scalar(@segs > 3)){
			$entry = BedCoord->new('chr' => $segs[0], 'start' => $segs[1], 'end' => $segs[2], 'other' => [@segs[3 .. $segNum]]);
		}else{
			$entry = BedCoord->new('chr' => $segs[0], 'start' => $segs[1], 'end' => $segs[2]);
		}
		my $bin = $binner->getbin($segs[1], $segs[2]);
		if(!($self->has_bed())){
			my %h;
			push(@{$h{$segs[0]}->{$bin}}, $entry);
			$self->bed(\%h);
		}else{
			my %h = %{$self->bed()};
			push(@{$h{$segs[0]}->{$bin}}, $entry);
			$self->bed(\%h);
		}
		$self->lines($self->lines() + 1);
	}
	close IN;
	print STDERR "[perlBed] Finished loading bed file: $file\n";
}

# Assumes the bed coordinates are merged
# Returns an integer value of the summed length of all bedcoords in this container
sub calculateBedLength{
	my ($self, $chr) = @_;
	
	my $sum = 0;
	foreach my $b (keys %{$self->bed()->{$chr}}){
		foreach my $bed (@{$self->bed()->{$chr}->{$b}}){
			$sum += $bed->end - $bed->start;
		}
	}
	return $sum;
}

# Returns BedCoord object if it intersects
# Returns '0' if no intersection
sub firstIntersect{
	my ($self, $chr, $start, $end) = @_;
	
	my $binner = kentBinTools->new();
        unless($self->has_bed()){
                print STDERR "[perlBed] Error! Called data_intersects before loading data!\n";
                exit;
        }
        if(!(exists($self->bed()->{$chr}))){
                return 0;
        }

        my @bins = $binner->searchbins($start, $end);
        foreach my $b (@bins){
                if(!(exists($self->bed()->{$chr}->{$b}))){
                        next;
                }
                foreach my $gc (@{$self->bed()->{$chr}->{$b}}){
                        if($binner->overlap($gc->start(), $gc->end(), $start, $end) > 1){
                                return $gc;
                        }
                }
        }
	return 0;
}

# Bool: does the input information intersect with the stored data
sub intersects {
	my ($self, $chr, $start, $end) = @_;
	
	my $binner = kentBinTools->new();
	unless($self->has_bed()){
		print STDERR "[perlBed] Error! Called data_intersects before loading data!\n";
		exit;
	}
	if(!(exists($self->bed()->{$chr}))){
		return 0;
	}
	
	my @bins = $binner->searchbins($start, $end);
	foreach my $b (@bins){
		if(!(exists($self->bed()->{$chr}->{$b}))){
			next;
		}
		foreach my $gc (@{$self->bed()->{$chr}->{$b}}){
			if($binner->overlap($gc->start(), $gc->end(), $start, $end) > 1){
				return 1;
			}
		}
	}
	return 0;
}

__PACKAGE__->meta->make_immutable;

1;
