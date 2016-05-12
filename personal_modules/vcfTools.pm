#!/usr/bin/perl
# A list of utilities designed to segment and split VCF files
# Right now, assumes that the file is indexed by bcftools -- will add functionality later in order to assert this fact

package VCFFilterAndSplit;
use Mouse;
use namespace::autoclean;

has ['vcffile', 'output'] => (is => 'ro', isa => 'Str', required => 1);
has 'ucsc' => (is => 'rw', isa => 'Str', predicate => 'has_region');
has 'indivs' => (is => 'rw', isa => 'Str', predicate => 'has_indivs');

sub generateSummaryTabFile{
	my ($self) = @_;
	my $file = $self->vcffile;
	my $output = $self->output;
	
	# Determine which individuals should be saved
	my @chosen = $self->_determineIndivCols;
	
	# Determine if a region should be filtered
	my $ucsc = ($self->has_region)? $self->ucsc : "";

	open(my $IN, "bcftools view $file $ucsc |") || die "Could not open bcftools output in main routine!\n";
	open(my $OUT, "> $output");

	while(my $line = <$IN>){
		if($line =~ /^##/){next;}
		elsif($line =~/^#CHROM/){
			# Identified the header line, now to generate the header for the output file
			chomp $line;
			my @segs = split(/\s+/, $line);
			print {$OUT} "CHR\t$segs[1]\t$segs[3]\t$segs[4]\t$segs[5]\tTYPE\tMUTATION\tPRIORITY\tGENE\tAA\t";
			my @data = $self->_returnDeterminedSegs(\@segs, \@chosen);
			print {$OUT} join("\t", @data) . "\n";
		}else{
			# A regular vcf line to process
			chomp $line;
			my @segs = split(/\t/, $line);
			print {$OUT} "$segs[0]\t$segs[1]\t$segs[3]\t$segs[4]\t$segs[5]\t";
			
			my @infosegs = split(/;/, $segs[7]);
			if($infosegs[0] eq "INDEL" || length($segs[3]) > 1 || length($segs[4]) > 1){
				print {$OUT} "INDEL";
			}else{
				print {$OUT} "SNP";
			}
			
			my @mutation;
			my @priority;
			my @gene;
			my @aas;
			my ($m, $p, $g, $aa);
			foreach my $i (@infosegs){
				if($i =~ /ANN=(.+$)/){
					my @effsegs = split(/\|/, $1);
					push(@mutation, $effsegs[1]);
					push(@aas, $effsegs[10]);
					push(@priority, $effsegs[2]);
					push(@gene, "$effsegs[3];$effsegs[4]");
				}
			}


			print {$OUT} "\t" . join(';', @mutation) . "\t" . join(';', @priority) . "\t" . join(';', @gene) . "\t" . join(';', @aas);
			
			# Get genotype entries and segregate them to just the calls
			my @data = $self->_returnDeterminedSegs(\@segs, \@chosen);
			foreach my $col (@data){
				my @forsegs = split(/:/, $col);
				#my @gsegs = split(/[\/\|]/, $forsegs[0]);

				print {$OUT} "\t$forsegs[0]";
			}
			print {$OUT} "\n";
		}
	}
	close $IN;
	close $OUT;
}

sub _returnDeterminedSegs{
	my ($self, $segRef, $saveRef) = @_;
	# Routine to slice the reference if needed
	my @data;
	foreach my $num (@{$saveRef}){
		push(@data, $segRef->[$num]);
	}
	return @data;
}

sub _determineIndivCols{
	my ($self) = @_;
	my $file = $self->vcffile;
	my @save;
	
	open(my $IN, "bcftools view $file |");
	my @chrline;
	while(my $line = <$IN>){
		if($line =~ /^##/){
			next;
		}elsif($line =~ /^#CHROM/){
			@chrline = split("\t", $line);
			last;
		}else{
			print STDERR "Error with malformed VCF file! Expected header but did not find one! Exiting...\n";
			exit;
		}
	}
	close $IN;
	
	# Determine if we need to slice the columns
	if($self->has_indivs){
		my $indivs = $self->indivs;
		my %indhash;
		open(my $IN, "< $indivs") || die "Could not open individual list file!\n";
		while(my $line = <$IN>){
			chomp $line;
			$indhash{$line} = 1;
		}
		close $IN;
		
		for(my $x = 9; $x < scalar(@chrline); $x++){
			if(exists($indhash{$chrline[$x]})){
				push(@save, $x);
			}
		}
	}else{
		# Just grab everything
		for(my $x = 9; $x < scalar(@chrline); $x++){
			push(@save, $x);
		}
	}
	
	return @save;
}

__PACKAGE__->meta->make_immutable;

1;