#!/usr/bin/perl
# This script takes an annotated bcf file (with index) and converts it into a tab delimited format for easy referencing
# Added vcf option to bypass bcf processing and allowed output of all variant lines

use strict;
use Getopt::Std;

my $usage = "perl $0 (-b <indexed bcf file> || -v <vcf file>) (-a <all> || -g <grep string> || -r <region in ucsc format>) -o <output file>\n";
my %opts;
getopt('brogv', \%opts);

unless((defined($opts{'b'}) || defined($opts{'v'})) && (defined($opts{'a'}) || defined($opts{'g'}) || defined($opts{'r'})) && defined($opts{'o'})){
	print $usage;
	exit;
}

if(defined($opts{'v'}) && defined($opts{'r'})){
	print "Error! Select indexed bcf for this option!\n$usage";
	exit;
}

open(OUT, "> $opts{o}");

my $vcf = defined($opts{'b'})? $opts{'b'} : $opts{'v'};

# Retrieve header sequence
my @header = getHeaderInfo($vcf);
print OUT join("\t", @header) . "\n";

#CHROM  POS     ID      REF     ALT    QUAL     FILTER  INFO    FORMAT  BIBR02.1        BIBR03.1        
my $in;
my $cmd = defined($opts{'b'})? "bcftools view $opts{b}" : "cat $opts{v}";
if(defined($opts{'r'})){
	open($in, "bcftools view $opts{b} $opts{r} |") || die "Could not open bcftools output!\n $usage";
}elsif(defined($opts{'g'})){
	open($in, "$cmd | grep \'$opts{g}\' |") || die "Could not open bcftools output!\n $usage";
}else{
	open($in, "$cmd |") || die "Could not open file output!\n $usage";
}

#my @header;
while(my $line = <$in>){
	if($line =~ /^#/){next;}
	else{
		my @saved = processVCFLine($line);
		print OUT join("\t", @saved) . "\n";
	}
}
close $in;
close OUT;

exit;

sub processVCFLine{
	my ($line) = @_;
	chomp $line;
	my @saved;
	
	my @segs = split(/\t/, $line);
	#print OUT "$segs[0]\t$segs[1]\t$segs[3]\t$segs[4]\t$segs[5]\t";
	push(@saved, ($segs[0], $segs[1], $segs[3], $segs[4], $segs[5]));
	#INDEL;IS=4,1.000000;DP=84;VDB=1.008941e-02;AF1=0.9542;AC1=138;DP4=2,2,52,8;MQ=28;FQ=-8.94;PV4=0.11,1,0.19,1;EFF=INTERGENIC(MODIFIER||||||||||1)
	my @infosegs = split(/;/, $segs[7]);
	if($infosegs[0] eq "INDEL"){
		#print OUT "INDEL";
		push(@saved, "INDEL");
	}else{
		#print OUT "SNP";
		push(@saved, "SNP");
	}

	#EFF=SYNONYMOUS_CODING(LOW|SILENT|acG/acA|T334|581|PRLR|protein_coding|CODING|ENSBTAT00000014437|9|1)
	#ANN=A|stop_gained|HIGH|TBC1D4|ENSBTAG00000005760|transcript|ENSBTAT00000066234|protein_coding|20/23|c.3217A>T|p.Arg1073*|3217/4499|3217/3915|1073/1304||,A|stop_gained|HIGH|TBC1D4|ENSBTAG00000005760|transcript|ENSBTAT00000007574|protein_coding|18/21|c.3211A>T|p.Arg1071*|3211/4493|3211/3909|1071/1302||,A|stop_gained|HIGH|TBC1D4|ENSBTAG00000005760|transcript|ENSBTAT00000065132|protein_coding|15/17|c.2473A>T|p.Arg825*|2473/2913|2473/2748|825/915||WARNING_TRANSCRIPT_NO_STOP_CODON;LOF=(TBC1D4|ENSBTAG00000005760|3|1.00)
	my $effindex = 0;
	for(my $x = 0; $x < scalar(@infosegs); $x++){
		if($infosegs[$x] =~ /^ANN=/){
			$effindex = $x; last;
		}
	}
	my @annosegs = split(/,/, $infosegs[$effindex]);
	my @mutation;
	my @priority;
	my @gene;
	my @aas;
	my ($m, $p, $g, $aa);
	foreach my $a (@annosegs){
		$a =~ s/ANN=//g;
		#($m) = $a =~ /(.+)\(.*\)/;
		my ($str) = $a =~ /.+\((.*)\)/;
		my @effsegs = split(/\|/, $str);
		$m = $effsegs[1];
		$p = $effsegs[2];
		$g = $effsegs[3];
		$aa = $effsegs[10];
		push(@aas, $aa);
		push(@mutation, $m);
		push(@priority, $p);
		push(@gene, $g);
	}

	#print OUT "\t" . join(';', @mutation) . "\t" . join(';', @priority) . "\t" . join(';', @gene) . "\t" . join(';', @aas);

	push(@saved, (join(';', @mutation), join(';', @priority), join(';', @gene), join(';', @aas)));
	for(my $x = 9; $x < scalar(@segs); $x++){
		my @forsegs = split(/:/, $segs[$x]);
		my @gsegs = split(/[\/\|]/, $forsegs[0]);

		#print OUT "\t" . join(';', @gsegs);
		push(@saved, join(';', @gsegs));
	}
	#print OUT "\n";
	
	return @saved;
}

sub getHeaderInfo{
	my ($file) = @_;
	my $IN; 
	my $cmd = ($file =~ /.+.vcf$/)? "cat $file | grep \'#CHROM\' |" : "bcftools view $file | grep \'#CHROM\' |";
	open($IN, "$cmd");
	my @header;
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\s+/, $line);
		push(@header, ("CHR", $segs[1], $segs[3], $segs[4], $segs[5], "TYPE", "MUTATION", "PRIORITY", "GENE", "AA"));
		push(@header, @segs[9 .. $#segs]);
		
		#print OUT "CHR\t$segs[1]\t$segs[3]\t$segs[4]\t$segs[5]\tTYPE\tMUTATION\tPRIORITY\tGENE\tAA\t";
		#print OUT join("\t", @segs[9 .. $#segs]) . "\n";
	}
	close $IN;
	return @header;
}