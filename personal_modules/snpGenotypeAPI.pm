#!/usr/bin/perl
# This is a collection of my classes and structures for use in genotype deconvolution

package PedMapPrinter;
use Mouse;
use namespace::autoclean;
use strict;

sub printOutPedMap{
	my($self, $genotypearray, $baseoutname, $genokey) = @_;
	open(MAP, "> $baseoutname.map");
	open(PED, "> $baseoutname.ped");
	my $fh;
	if(defined($genokey)){
		open($fh, "< $genokey") || die "[PRINTER] Could not open genokey file: $genokey!\n";
	}
	
	my $mapprinted = 0;
	foreach my $genotype (@{$genotypearray}){
		if(ref($genotype) eq "AIPLGenotype"){
			my %aipl = ('0' => 'B B', '1' => 'A B', '2' => 'A A', '3' => '0 0', '4' => '0 0', '5' => '0 0');
			if(defined($genokey)){
				my $code = <$fh>;
				my @codes = split(//, $code);
				$aipl{'0'} = "$codes[1] $codes[1]";
				$aipl{'1'} = "$codes[0] $codes[1]";
				$aipl{'2'} = "$codes[0] $codes[0]";
			}
			#print MAP "$cn\t$segs[0]_cnregion_$cnvritr\t$segs[1]\t$segs[2]\n";
			#unshift(@matrix, ($breeds{$sub}, $an, 0, 0, 1, $pops{$sub}));
			my @chrs = $genotype->maplocs->get_chrs();
			my @locs = $genotype->maplocs->get_locs();
			my @genos = $genotype->geno->get_array();
			
			if(!$mapprinted){
				for(my $x = 0; $x < scalar(@chrs); $x++){
					my $c = $chrs[$x];
					my $l = $locs[$x];
					my ($cn) = $c =~ /chr(.+)/;
					print MAP "$cn\t$x\_loc\t0\t$l\n";
				}
			}
			
			print PED $genotype->population() . " " . $genotype->anname() . " 0 0 1 0";
			foreach my $g (@genos){
				my $conv = $aipl{$g};
				print PED " $conv";
			}
			print PED "\n";
		}else{
			die "[PEDMAP PRINTER] Cannot handle class: " . ($genotype->blessed()) . " currently!\n";
		}
		$mapprinted = 1;
	}
	close PED;
	close MAP;
	if(defined($genokey)){
		close $fh;
	}
}

__PACKAGE__->meta->make_immutable;

package CreateGenotypeCode;
use Mouse;
use namespace::autoclean;

sub createGenoCodeMap{
	my ($self, $file, $mapbed, $genofile) = @_;
	my $tfile1 = "mylocs_" . (int(rand(1000000))) . ".bed";
	$mapbed->createBedFile($tfile1);
	
	open(OUT, "> $genofile");
	open(IN, "intersectBed -a $file -b $tfile1 | uniq |") || die "[CREATECODE] Could not access bedtools to create code file!\n";
	# to get rid of non-unique lines
	my ($prechr, $prestart, $preend);
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($prechr eq $segs[0] && $prestart == $segs[1] && $preend == $segs[2]){next;}
		print OUT "$segs[3]\n";
		$prechr = $segs[0]; $prestart = $segs[1]; $preend = $segs[2];
	}
	close IN;
	close OUT;
	system("rm $tfile1");
}

__PACKAGE__->meta->make_immutable;

package AIPLGenotype;
use Mouse;
use namespace::autoclean;
use strict;

has ['anname', 'population'] => (is => 'rw', isa => 'Str');
has ['ankey', 'sex', 'snpnum'] => (is => 'rw', isa => 'Int');

has 'geno' => (is => 'rw', isa => 'genotype');
has 'maplocs' =>(is => 'rw', isa => 'mapBed');

sub convertGenoStr{
	my ($self, $gstr, $population) = @_;
	my @segs = split(/\s+/, $gstr);
	my @gsegs = split(//, $segs[4]);
	
	$self->anname($segs[0]);
	$self->ankey($segs[1]);
	$self->sex($segs[2]);
	$self->snpnum($segs[3]);
	
	my $genotype = genotype->new();
	$genotype->genotype(\@gsegs);
	
	$self->geno($genotype);
	$self->population($population);
}



__PACKAGE__->meta->make_immutable;

package mapBed;
use Mouse;
use namespace::autoclean;
use strict;

has 'chrarray' => (
	traits => ['Array'],
	is => 'rw',
	isa => 'ArrayRef[Str]',
	default => sub{[]},
	handles =>{
		'push_chr' => 'push',
		'get_chrs' => 'elements',
		'getc' => 'get',
		'countc' => 'count',
	}
);

has 'locarray' => (
	traits => ['Array'],
	is => 'rw',
	isa => 'ArrayRef[Int]',
	default => sub{[]},
	handles =>{
		'push_loc' => 'push',
		'get_locs' => 'elements',
		'getl' => 'get',
		'countl' => 'count',
	}
);

sub convertFromBed{
	my ($self, $file) = @_;
	open(IN, "< $file") || die "[MAPBED - BED] Could not open input bed file!\n";
	my (@chrs, @locs);
	while(my $line = <IN>){
		chomp $line;
		$line =~ s/\r//g;
		my @segs = split(/\t/, $line);
		push(@chrs, $segs[0]);
		push(@locs, $segs[1]);
		#$self->push_chr($segs[0]);
		#$self->push_loc($segs[1]);
	}
	$self->chrarray(\@chrs);
	$self->locarray(\@locs);
	print STDERR "[MAPBED - BED] Finished loading!\n";
	close IN;
}

sub convertFromCHRDATA{
	my ($self, $file) = @_;
	open(IN, "< $file") || die "[MAPBED - CHRDATA] Could not open input chromosome.data file!\n";
	my $h = <IN>;
	my (@chrs, @locs);
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\s+/, $line);
		my $chr = "chr" . $segs[1];
		push(@chrs, $chr);
		push(@locs, $segs[4]);
		#$self->push_chr($chr);
		#$self->push_loc($segs[4]);
	}
	$self->chrarray(\@chrs);
	$self->locarray(\@locs);
	print STDERR "[MAPBED - CHRDATA] Finished loading!\n";
	close IN;
}

sub removeLocs{
	my ($self, $mapbed, $genotype) = @_;
	my $tfile1 = "mylocs_" . (int(rand(1000000))) . ".bed";
	my $tfile2 = "yourlocs_" . (int(rand(1000000))) . ".bed";
	$self->createBedFile($tfile1, $genotype);
	$mapbed->createBedFile($tfile2);
	
	my (@newchrs, @newlocs, @newgenos);
	my $newGeno = genotype->new();
	my $newMap = mapBed->new();
	# to get rid of non-unique lines
	my ($prechr, $prestart, $preend);
	open(IN, "intersectBed -a $tfile1 -b $tfile2 -v |") || die "[MAPBED - RMLOCS] Could not access intersectBed output for $tfile1 and $tfile2!\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($prechr eq $segs[0] && $prestart == $segs[1] && $preend == $segs[2]){next;}
		push(@newchrs, $segs[0]);
		push(@newlocs, $segs[1]);
		if(defined($genotype)){
			#$newGeno->push_geno($segs[3]);
			push(@newgenos, $segs[3]);
		}
		$prechr = $segs[0]; $prestart = $segs[1]; $preend = $segs[2];
	}
	close IN;
	system("rm $tfile1 $tfile2");
	
	#$self->locarray(\@newlocs);
	#$self->chrarray(\@newchrs);
	$newMap->locarray(\@newlocs);
	$newMap->chrarray(\@newchrs);
	if(defined($genotype)){
		$newGeno->genotype(\@newgenos);
		return $newMap, $newGeno;
	}
	return $newMap;
}

sub intersectLocs{
	my ($self, $mapbed, $genotype) = @_;
	my $tfile1 = "mylocs_" . (int(rand(1000000))) . ".bed";
	my $tfile2 = "yourlocs_" . (int(rand(1000000))) . ".bed";
	$self->createBedFile($tfile1, $genotype);
	$mapbed->createBedFile($tfile2);
	
	my (@newchrs, @newlocs, @newgenos);
	my $newGeno = genotype->new();
	my $newMap = mapBed->new();
	# to get rid of non-unique lines
	my ($prechr, $prestart, $preend);
	open(IN, "intersectBed -a $tfile1 -b $tfile2 |") || die "[MAPBED - INLOCS] Could not access intersectBed output for $tfile1 and $tfile2!\n";
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($prechr eq $segs[0] && $prestart == $segs[1] && $preend == $segs[2]){next;}
		push(@newchrs, $segs[0]);
		push(@newlocs, $segs[1]);
		if(defined($genotype)){
			#$newGeno->push_geno($segs[3]);
			push(@newgenos, $segs[3]);
		}
		$prechr = $segs[0]; $prestart = $segs[1]; $preend = $segs[2];
	}
	close IN;
	system("rm $tfile1 $tfile2");
	
	#$self->locarray(\@newlocs);
	#$self->chrarray(\@newchrs);
	$newMap->locarray(\@newlocs);
	$newMap->chrarray(\@newchrs);
	if(defined($genotype)){
		$newGeno->genotype(\@newgenos);
		return $newMap, $newGeno;
	}
	return $newMap;
}

sub createBedFile{
	my ($self, $filename, $genotype) = @_;
	open(OUT, "> $filename") || die "[MAPBED - CRBEDFILE] Could not create output bed file: $filename!\n";
	
	if(defined($genotype)){
		if($self->countl() != $genotype->countg()){
			die "[MAPBED - CRBEDFILE] Count of map locations does not match genotype count for output file: $filename!\n";
		}
		
		for(my $x = 0; $x < $self->countl(); $x++){
			my $chr = $self->getc($x);
			my $start = $self->getl($x);
			my $end = $start + 1;
			my $gen = $genotype->getg($x);
			print OUT "$chr\t$start\t$end\t$gen\n";
		}
	}else{
		for(my $x = 0; $x < $self->countl(); $x++){
			my $chr = $self->getc($x);
			my $start = $self->getl($x);
			my $end = $start + 1;
			my $gen = $x;
			print OUT "$chr\t$start\t$end\t$gen\n";
		}
	}
	close OUT;
}

__PACKAGE__->meta->make_immutable;

package genotype;
use Mouse;
use namespace::autoclean;
use strict;

has 'genotype' => (
	traits => ['Array'],
	is => 'rw',
	isa => 'ArrayRef[Int]',
	default => sub{[]},
	handles => {
		'push_geno' => 'push',
		'get_array' => 'elements',
		'getg' => 'get',
		'countg' => 'count',
	},
);

__PACKAGE__->meta->make_immutable;

1;