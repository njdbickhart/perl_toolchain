#!/usr/bin/perl
# This script is a one-off script designed to process a series of scaffolds against an ordered RH map to generate a larger AGP file
# My strategy is to create a tree map of contigs from the agp file and then associate them, in order, with the RH map

use strict;
use Getopt::Std;

my %opts;

my $usage = "perl $0 -a inputAGPFile -r inputRHMapFile -o outputAGP\n";

getopt('aro', \%opts);

unless(defined($opts{'a'}) && defined($opts{'r'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

# Create the tree
my $worker = AGPTree->new();
$worker->addAGP($opts{'a'});

# Process the data
$worker->processRHmap($opts{'r'}, $opts{'o'});

exit;

BEGIN{
package AGPNode;
use Mouse;
use namespace::autoclean;

has ['scafstart', 'scafend', 'qstart', 'qend', 'ovpstart', 'ovpend'] => (is => 'rw', isa => 'Int');
has ['qname', 'orient', 'scaf', 'chr'] => (is => 'rw', isa => 'Str');

sub processLine{
	my ($self, $line) = @_;
	chomp($line);
	$line =~ s/\r//g;
	
	my @segs = split(/\t/, $line);
	$self->scafstart($segs[1]);
	$self->scafend($segs[2]);
	$self->qstart($segs[6]);
	$self->qend($segs[7]);
	#$self->ovpstart($segs[9]);
	#$self->ovpend($segs[10]);
	
	my @qsegs = split(/\_/, $segs[5]);
	$self->scaf($segs[0]);
	$self->qname($qsegs[0]);
	$self->orient($segs[8]);
}

sub getArray{
	my ($self) = @_;
	my @data = ($self->chr, $self->scafstart, $self->scafend, 0, "W", $self->qname, $self->qstart, $self->qend, $self->orient, $self->scaf);
	return @data;
}

__PACKAGE__->meta->make_immutable;

package AGPTree;
use Mouse;
use namespace::autoclean;

# Tree hashkeys are the superscaffold IDs
has 'tree' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub {{}},
	handles => {
		'addnode' => 'set',
		'getref' => 'get',
	});

# The Lookup table accepts a contig ID and returns the appropriate superscaffold ID
has 'lookup' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub {{}},
	handles => {
		'addcontig' => 'set',
		'getcontig' => 'get',
		'cexists' => 'exists',
	});
	
sub addAGP{
	my ($self, $file) = @_;
	open(IN, "< $file") || die "Could not open agp file!\n";
	my %holder;
	while(my $line = <IN>){
		if($line =~ /^#/){next;}
		chomp $line;
		my @segs = split(/\t/, $line);
		my @qsegs = split(/\_/, $segs[5]);
		if($segs[4] ne "W"){next;}
		my $node = AGPNode->new();
		$node->processLine($line);
		$self->addcontig($qsegs[0] => $segs[0]);
		push(@{$holder{$segs[0]}}, $node);
	}
	close IN;
	$self->tree(\%holder);
}

sub processRHmap{
	my ($self, $rh, $out) = @_;
	open(IN, "< $rh") || die "Could not open RH map file!\n";
	open(OUT, "> $out");
	
	my %used;
	my $h = <IN>;
	my @holder;
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		
		my @qsegs = split(/\_/, $segs[5]);
		
		if($segs[5] =~ /\*/ || $segs[5] =~ /\#N\/A/){next;}
		
		$segs[5] = $qsegs[0];
		push(@holder, \@segs);
	}
	close IN;
	
	my @working;
	my $prev_scaf = [];
	my $curchr = "1"; #Assuming we always start with chr 1
	for(my $x = 0; $x < scalar(@holder); $x++){
		my $ctg = $holder[$x]->[5];
		if($holder[$x]->[5] eq '*' || $holder[$x]->[5] eq '#N/A' || $holder[$x]->[5] eq ''){
			next;
		}
		
		if($x == 2289){
			print "start!\n";
		}
		
		if($curchr ne $holder[$x]->[1] && scalar(@working) > 0 && scalar(@{$prev_scaf}) > 0){
			# Chromosome ended -- time to work on the data we have
			$curchr = $holder[$x]->[1];
			my @lines = $self->processArrays($prev_scaf, \@working);
			foreach my $l (@lines){
				print OUT join("\t", @{$l});
				print OUT "\n";
			}
			@working = ();
		}elsif($self->cexists($ctg)){
			# If the contig is in a super scaffold, find all associated contigs
			
			my %ctgSet;
			my $scaf = $self->getcontig($ctg);
			my $scaf_array = $self->getref($scaf);
			$prev_scaf = $scaf_array;
			foreach my $s (@{$scaf_array}){
				$ctgSet{$s->qname} = 1;
			}
			
			if(exists($ctgSet{"utg60833"})){
				print "start!\n";
			}
			
			for(my $y = $x; $y < scalar(@holder); $y++){
				if($y == 2289){
					print "hey!\n";
				}
				if(exists($ctgSet{$holder[$y]->[5]})){
					push(@working, $holder[$y]);
				}else{
					# Race condition to try to find if a contig exists in the middle of two known scaffolded contigs
					my $sandwich = 0;
					my $probctg = $holder[$y]->[5];
					for(my $z = $y + 1; $z < scalar(@holder); $z++){
						if($holder[$z]->[5] ne $probctg && exists($ctgSet{$holder[$z]->[5]})){
							# contig is sandwhiched between two parts of the super scaffolds
							push(@working, $holder[$y]);
							$sandwich = 1;
							last;
						}elsif($holder[$z]->[5] ne $probctg && !(exists($ctgSet{$holder[$z]->[5]}))){
							# There are two contigs that are not part of the super scaffold! Abort!
							last;
						}
					}
					if(!$sandwich){
						# We reached the end!
						last;
					}
				}
				$x = $y;
			}
			
			# Now to process the arrays
			my @lines = $self->processArrays($scaf_array, \@working);
			foreach my $l (@lines){
				print OUT join("\t", @{$l});
				print OUT "\n";
			}
			@working = ();
		}else{
			
			# This was a contig not joined by the BioNano scaffolding -- printing out in a different format for the AGP
			my $cur = $holder[$x];
			#my @values = ($cur->[1], $cur->[2], $cur->[2], $x, "W", $cur->[5], $cur->[6], $cur->[6], "*");
			#print OUT join("\t", @values);
			#print OUT "\n";
		}
	}
	close OUT;

}

sub processArrays{
	my ($self, $agparray, $rharray) = @_;
	# OK, since the AGP has orientation, and the RH map does not, I need to align the AGP in the right orientation to print it out
	
	# account for simple alignment of one contig
	my $aidsame = ($agparray->[0]->qname eq $agparray->[-1]->qname)? 1 : 0;
	my $ridsame = ($rharray->[0]->[5] eq $rharray->[-1]->[5])? 1 : 0;
	
	# Set all of the agp chromosomes to the RH chr
	my $chr = $rharray->[0]->[1];
	foreach my $agp (@{$agparray}){
		$agp->chr($chr);
	}
	my %ctgSet;
	foreach my $s (@{$agparray}){
		$ctgSet{$s->qname} = 1;
	}
	if(exists($ctgSet{"utg60833"})){
		print "start!\n";
	}
	my @container;
	if($aidsame && $ridsame){
		if(scalar(@{$rharray}) > 1){
			# there are more than one points of coverage for the RH map -- we can orient the direction of the agp file
			my ($rorientrev, $aorientrev) = $self->determineOrient($rharray->[0]->[6], $rharray->[1]->[6], $agparray->[0]->orient);
			if(($rorientrev && $aorientrev) || (!$rorientrev && !$aorientrev)){
				# Same orientation
				foreach my $agp (@{$agparray}){
					my @data = $agp->getArray;
					push(@container, \@data);
				}
			}else{
				# I'm going to flip the agp orientation by default here
				foreach my $agp (@{$agparray}){
					$agp->orient(($agp->orient eq "-")? "+" : "-");
					my @data = $agp->getArray;
					push(@container, \@data);
				}
			}
			return @container;
		}
		# there was only one entry
		foreach my $agp (@{$agparray}){
			my @data = $agp->getArray;
			push(@container, \@data);
		}
		return @container;
	}
	
	# Now we've gotta get more complicated
	
	my $revcontig = ($rharray->[0]->[5] eq $agparray->[-1]->qname && $rharray->[-1]->[5] eq $agparray->[0]->qname)? 1 : 0;
	
	# Run through the rh array to find entries that may have been missing from the agp file and add them in
	my @ctgorder;
	foreach my $r (@{$rharray}){
		if(scalar(@ctgorder) == 0 && $r->[5] ne ''){
			push(@ctgorder, $r);
		}else{
			if($ctgorder[-1]->[5] ne $r->[5] && $r->[5] ne ''){
				push(@ctgorder, $r);
			}
		}
	}
	
	my $ctgorderidx = ($revcontig)? scalar(@ctgorder) - 1 : 0;
	for(my $x = ($revcontig)? scalar(@{$agparray}) - 1 : 0; ;){
		if($revcontig && $x == 0){
			last;
		}elsif(!$revcontig && $x >= scalar(@{$agparray})){
			last;
		}
		
		my $curagp = $agparray->[$x];
		
		# Check if this agp map contig is preceeded by a different rh map contig
		# There should be only one but I will check
		#if($ctgorder[$ctgorderidx]->[5] ne $curagp->qname){
		#	my $cur = $ctgorder[$ctgorderidx];
		#	my @values = ($cur->[1], $cur->[2], $cur->[2], $x, "W", $cur->[5], $cur->[6], $cur->[6], "*");
		#	push(@container, \@values); 
		#	if($revcontig){
		#		$ctgorderidx--;
		#	}else{
		#		$ctgorderidx++;
		#	}
		#	if($ctgorder[$ctgorderidx]->[5] ne $curagp->qname){
		#		print STDERR "ISSUE with contig order! " . $ctgorder[$ctgorderidx]->[5] . "\tvs\t" . $curagp->qname . "\n";
		#	}
		#}
		
		if($revcontig){
			# Flip the sign because the ordering is different
			$curagp->orient(($curagp->orient eq "+")? "-" : "+");
		}
		
		my @data = $curagp->getArray;
		push(@container, \@data);
		
		if($revcontig){
			$x--;
		}else{
			$x++;
		}
	}
	
	return @container;
}

sub determineOrient{
	my ($self, $rhcoord1, $rhcoord2, $aorient) = @_;
	my $rorientrev = ($rhcoord1 > $rhcoord2)? 1 : 0;
	my $aorientrev = ($aorient eq "-")? 1 : 0;
	return ($rorientrev, $aorientrev);
}

1;
}