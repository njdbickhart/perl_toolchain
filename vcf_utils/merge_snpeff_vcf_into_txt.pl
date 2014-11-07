#!/usr/bin/perl
# This is a script designed to take two filtered snpeff vcfs and to combine all of the variants between them to generate one table

use strict;
use Class::Struct;

struct('snp' => {
	'pos' => '$',
	'type' => '@',
	'eff' => '@',
	'gene' => '$',
	'animals' => '@',
});

chomp(@ARGV);
my $usage = "usage: $0 <filtered vcf> ...\n";

if(scalar(@ARGV) < 1){
	print $usage;
	exit;
}

my %snpholder;
foreach my $file (@ARGV){
	my $name = '';
	open(IN, "< $file") || die "problem opening file: $file!\n$usage";
	while(my $line = <IN>){
		chomp $line;
		if($line =~ /^##/){next;}
		elsif($line =~ /^#CHROM/){
			my @segs = split(/\s+/, $line);
			$name = $segs[9];
		}else{
			if($name eq ''){ die "could not get name of animal!\n$usage";}
			my @segs = split(/\s+/, $line);
			my ($eff) = $segs[7] =~ m/;EFF=(.{3,25})\(/;
			my ($type) = $segs[7] =~ m/;TYPE=(.{3,5});/;
			my ($modstr) = $segs[7] =~ m/\(MODIFIER(|.*|.*|.*|.*|.*|.*|.*|.*|)\)/;
			my @modsegs = split(/\|/, $modstr);
			my $gene = $modsegs[5];
			
			my @gsegs = split(/:/, $segs[9]);
			my $allele;
			if($gsegs[0] eq "0/1"){
				$allele = "het";
			}elsif($gsegs[0] eq "1/1"){
				$allele = "hom";
			}else{
				$allele = "NA";
			}
			
			my $pos = $segs[1];
			if(exists($snpholder{$pos})){
				push(@{$snpholder{$pos}->type()}, $type);
				push(@{$snpholder{$pos}->eff()}, $eff);
				push(@{$snpholder{$pos}->animals()}, "$name:$allele");
			}else{
				$snpholder{$pos} = snp->new('pos' => $pos, 'type' => [$type], 'eff' => [$eff], 'gene' => $gene, 'animals' => ["$name:$allele"]);
			}
		}
	}
	close IN;

}

print "chr\tposition\ttype\teffect\tgene\tanimals\t#animals\n";
foreach  my $p (sort{$a <=> $b} keys(%snpholder)){
	my $g = $snpholder{$p};
	print "chr20\t" . $g->pos() . "\t" . join(";", @{$g->type()}) . "\t" . join(";", @{$g->eff()}) . "\t" . $g->gene() . "\t" . join(";", @{$g->animals()}) . "\t" . scalar(@{$g->animals()}) . "\n";
}

exit;