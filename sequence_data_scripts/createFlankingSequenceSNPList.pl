#!/usr/bin/perl
# This script is designed to take an input tab delimited SNP list, and generate a 150 bp flanking sequence file for probe design
# Probe list: chr, pos, Ref, Alt
# Output format: SNP_ID, Sequence, Chr, Pos, Minor Allele, Notes, VariantPhred, 5'MapQ, 3'MapQ 

use strict;
use Getopt::Std;
our $flanking = 150;

my $usage = "perl $0 -f <input fasta> -p <probe site list> -v <gzipped vcfs<comma separated>> -m <mapq assoc tab file> -o <output tab>\n";
my %opts;
getopt('fpvom', \%opts);

unless(defined($opts{'f'}) && defined($opts{'p'}) && defined($opts{'v'}) && defined($opts{'o'}) && defined($opts{'m'})){
	print $usage;
	exit;
}

my %initTargSites; # Initial list of pos coords for SNP sites chosen by Doro
my %TargetSites; # array of varSite objects for target sites; {chr}->[varSites]
my %OvlpSites; # array of overlapping sites per SNP id; {chr}->{pos}->[varSites]

# Start by reading in the list of probe sites
open(my $IN, "< $opts{p}") || die "Could not open probe site list file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	push(@{$initTargSites{$segs[0]}}, [$segs[0], $segs[1], $segs[2], $segs[3]]);
}
close $IN;

# Run through VCFs, collecting overlapping sites and generating objects for target sites
my @vcfs = split(/,/, $opts{v});
foreach my $v (@vcfs){
	open(my $IN, "gunzip -c $v |");
	while(my $line = <$IN>){
		if($line =~ /^#/){next;} # skip comment line
		
		chomp $line;
		my @segs = split(/\t/, $line);
		# check for target site status and/or overlap status
		foreach my $j (@{$initTargSites{$segs[0]}}){
			# target site definition
			if($j->[1] == $segs[1] && $j->[3] eq $segs[4]){
				my $var = new varSite();
				$var->ConvertFromVCF($line);
				push(@{$TargetSites{$segs[0]}}, $var);
			}
			
			# overlap site definition
			my $ostart = $j->[1] - ($flanking + 1);
			my $oend = $j->[1] + ($flanking);
			if($segs[1] <= $oend && $segs[1] >= $ostart){
				my $var = new varSite();
				$var->ConvertFromVCF($line);
				push(@{$OvlpSites{$j->[0]}->{$j->[1]}}, $var);
			}
		}
	}
	close $IN;
}

#Reconciling overlapping variants with their target sites
foreach my $chr (keys(%TargetSites)){
	foreach my $v (@{$TargetSites{$chr}}){
		my $pos = $v->pos();
		foreach my $o (@{$OvlpSites{$chr}->{$pos}}){
			$v->addVar($o);
		}
		
		# Extract flanking sequence from Fasta file for each targetSite
		my $start = $pos - ($flanking + 1);
		my $end = $pos + $flanking;
		if($start < 0){$start = 0;}
		open(my $IN, "samtools faidx $opts{f} $chr:$start-$end |");
		my $head = <$IN>;
		my $seq;
		while(my $line = <$IN>){
			chomp $line;
			$seq .= $line;
		}
		close $IN;
		$v->InitSeq($seq);
		
		$v->GenerateFormatSeq();
	}
}

# Associate mapq scores with existing varsites
open(my $IN, "< $opts{m}") || die "Could not open mapq assoc tab file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	foreach my $v (@{$TargetSites{$segs[0]}}){
		if($v->ID eq $segs[6]){
			$v->fivePrime($segs[4]);
			$v->threePrime($segs[5]);
		}
	}
}
close $IN;

# Print output and exit
open(my $OUT, "> $opts{o}");
foreach my $chr (sort{$a cmp $b} keys(%TargetSites)){
	foreach my $v (sort{$a->pos <=> $b->pos} @{$TargetSites{$chr}}){
		my $str = $v->FormatOutput();
		print {$OUT} "$str\n";
	}
}

exit;

BEGIN{
package varSite;
use Mouse;
use namespace::autoclean;

our %iupac = ("AG" => "R", "CT" => "Y", "GC" => "S", "AT" => "W", "GT" => "K", "AC" => "M");
our $flanking = 150;

has ['ID', 'InitSeq', 'FormatSeq', 'Chr', 'Ref', 'Alt'] => (is => 'rw', isa => 'Str');
has ['pos', 'Qual', 'MAF', 'fivePrime', 'threePrime'] => (is => 'rw', isa => 'Num', default => 0);
has 'OvlpVars' => (is => 'rw', isa => 'ArrayRef[varSite]', default => sub{[]},
	handles => {
		'addVar' => 'push',
	});

sub ConvertFromVCF{
	my ($self, $vcfline) = @_;
	chomp($vcfline);
	my @segs = split(/\t/, $vcfline);
	$self->ID($segs[2]); 
	$self->Chr($segs[0]); 
	$self->Ref($segs[3]);
	my @altsegs = split(/,/, $segs[4]); 
	$self->Alt($altsegs[0]); 
	$self->pos($segs[1]); 
	$self->Qual($segs[5]);
	
	my ($ac, $an) = $segs[7] =~ /.*AC=(.{1,4});.*AN=(.{1,4});.*/;
	my @acsegs = split(/,/, $ac);
	
	$self->MAF($acsegs[0] / $an);
}
	
sub GenerateFormatSeq{
	my ($self)  = @_;
	my $OvlpVars = $self->OvlpVars();
	my $initSeq = $self->InitSeq();
	if(scalar(@{$OvlpVars}) < 1){
		$self->FormatSeq($self->_BracketSeq($initSeq, $self->Ref(), $self->Alt()));
	}else{
		# Add IUPAC alterations and then generate formatted sequence
		foreach my $vars (sort{$a->pos <=> $b->pos} @{$self->OvlpVars}){
			my $strPos = ($vars->pos - $self->pos) + $flanking + 1;
			if($strPos < 0){
				print "Error with relative position calculation!\n";
			}
			my $curRef = substr($initSeq, $strPos, 1); 
			my $curAlt = $vars->Alt;
			my @nucs = ($curRef, $curAlt);
			@nucs = sort{$a cmp $b} @nucs;
			my $value = join("", @nucs);
			if(exists($iupac{$value})){
				substr($initSeq, $strPos, 1) = $iupac{$value};
			}else{
				substr($initSeq, $strPos, 1) = "N";
			}
		}
		$self->FormatSeq($self->_BracketSeq($initSeq, $self->Ref(), $self->Alt()));
	}
}

sub _BracketSeq{
	my ($self, $sequence, $ref, $alt) = @_;
	my $format = $sequence;
	substr($format, $flanking + 1, 1) = "[$ref\/$alt]";
	return $format;
}

sub FormatOutput{
	my ($self) = @_;
	
	# Output format: SNP_ID, Sequence, Chr, Pos, Minor Allele, Notes, VariantPhred, 5'MapQ, 3'MapQ 
	my @segs = ($self->ID, $self->FormatSeq, $self->Chr, $self->pos, $self->MAF, "", $self->Qual, $self->fivePrime, $self->threePrime);
	my $str = join("\t", @segs);
	return $str;
}

1;
}