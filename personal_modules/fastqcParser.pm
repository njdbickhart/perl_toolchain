#!/usr/bin/perl
# This module is designed to be an encapsulated interpreter of a fastqc_data.txt summary file
# It will process the file, summarize the statistics and, on demand, provide a suitable tab delimited summary
# The whole goal is to condense the information from fastqc across an entire sample into a spreadsheet-style row

# Required in constructor:
# file => fastq file, sample => sample_name, library => library_name, readNum => first_or_second_read, log => simpleLogger

# Pipeline:
# 	$object->runFastqc($fastqc_exe)
#	$object->parseStats()
#	$object->cleanUp( (1))  <- optional removal of folders and (zip file)
#	my @headers = $object->getHeaderArray()
#	my @values = $object->getOutArray()  #Output is a one dimensional array of values

package fastqcParser;
use Mouse;
use namespace::autoclean;
use simpleLogger;
use File::Basename;

our @seHeaders = ("Sample", "Library", "FastqFile", "ReadNum", "TotalReads", "FilteredReads", "ReadLength", "TotalGCPerc", "BPQualityPass",
	"ReadQSPass", "BPGCContentPass", "BPNContentPass", "OverRepPass", "KmerPass", "SeqDupPerc",
	"First25Qual", "Second25Qual", "Third25Qual", "Fourth25Qual", "#ReadsQS2-11", "#ReadsQS12-21", 
	"#ReadsQS22-31", "#ReadsQS32-40", "First25GC", "Second25GC", "Third25GC", "Fourth25GC", 
	"First25N", "Second25N", "Third25N", "Fourth25N" );

has ['pass', 'overrep', 'kmer'] => (is => 'rw', isa => 'Str', default => 'PASS');
has ['file', 'sample', 'library'] => (is => 'ro', isa => 'Str', required => 1);
has 'log' => (is => 'ro', isa => 'simpleLogger', required => 1);
has 'readNum' => (is => 'ro', isa => 'Num', required => 1);
has 'seqlendist' => (is => 'rw', isa => 'Str');
has 'seqdup' => (is => 'rw', isa => 'Num');

# Locations of fastqc folder and zip
has ['folder', 'zip'] => (is => 'rw', isa => 'Str');

# Containers for individual stat modules
has 'sampstats' => (is => 'rw', isa => 'sampStats');
has 'pbpquality' => (is => 'rw', isa => 'pbpQuality');
has 'pseqquality' => (is => 'rw', isa => 'pseqQuality');
has 'pbpgccontent' => (is => 'rw', isa => 'pbpGCContent');
has 'pbpncontent' => (is => 'rw', isa => 'pbpNContent');

sub runFastqc{
	my ($self, $fastqc) = @_;
	if(!defined($fastqc)){
		$self->log->Fatal("[FQCPARSE]", "Did not receive fastqc executable!");
	}
	
	$self->log->Info("[FQCPARSE]", "Beginning fastqc runtime on file: " . $self->file);
	
	# Determine file output 
	my ($filename, $dirs, $ext) = fileparse($self->file);
	#if($ext eq 'gz'){
		my @fsegs = split(/\./, $filename);
		$self->folder("$dirs/$fsegs[0]" . "_fastqc");
		$self->zip("$dirs/$fsegs[0]" . "_fastqc.zip");
	#}else{
	#	$self->folder("$dirs/$filename" . "_fastqc");
	#	$self->zip("$dirs/$filename" . "_fastqc.zip");
	#}
	
	system("$fastqc -q " . $self->file);
	$self->log->Info("[FQCPARSE]", "Finished fastqc runtime on file: " . $self->file);
}

sub parseStats{
	my ($self) = @_;
	
	$self->log->Info("[FQCPARSE]", "Parsing fastqc stats on file: " . $self->file);
	
	open(IN, "< " . $self->folder . "/fastqc_data.txt") || $self->log->Fatal("[FQCPARSE]", "Error accessing fastqc_data.txt in folder: " . $self->folder . " for file: " . $self->file . "!");
	my $store;
	my $inloop = 0;
	while(my $line = <IN>){
		if($line =~ /^>>/ && $line !~ />>END_MODULE/){
			$inloop = 1;
			$store .= $line;
		}elsif($line =~ />>END_MODULE/){
			$inloop = 0;
			
			# TODO: rework this section to call individual class parsers rather than a blanket sub
			my @linesegs = split(/\n/, $store);
			my ($name, $pass) = $linesegs[0] =~ m/>>(.+)\s+(.+)$/;
			$store = '';

			if($name eq "Per base sequence quality"){
				my $pbpquality = pbpQuality->new();
				$pbpquality->fillContainers(\@linesegs, $pass);
				$self->pbpquality($pbpquality);
			}elsif($name eq "Per sequence quality scores"){
				my $pseqquality = pseqQuality->new();
				$pseqquality->fillContainers(\@linesegs, $pass);
				$self->pseqquality($pseqquality);
			}elsif($name eq "Per base sequence content"){
				my $pbpgccontent = pbpGCContent->new();
				$pbpgccontent->fillContainers(\@linesegs, $pass);
				$self->pbpgccontent($pbpgccontent);
			}elsif($name eq "Per base GC content" || $name eq "Per sequence GC content"){
			}elsif($name eq "Per base N content"){
				my $pbpncontent = pbpNContent->new();
				$pbpncontent->fillContainers(\@linesegs, $pass);
				$self->pbpncontent($pbpncontent);
			}elsif($name eq "Sequence Length Distribution"){
				my @lengths;
				for(my $x = 1; $x < scalar(@linesegs); $x++){
					if($linesegs[$x] =~ /^#/){next;}
					my @msegs = split(/\s+/, $linesegs[$x]);
					push(@lengths, "$msegs[0]=$msegs[1]");
				}	
				$self->seqlendist(join(";", @lengths));
			}elsif($name eq "Sequence Duplication Levels"){
				foreach my $v (@linesegs){
					if($v =~ /#Total Duplicate Percentage\s+(\d+)/){
						$self->seqdup($1);
					}
				}
			}elsif($name eq "Overrepresented sequences"){
				$self->overrep($pass);
			}elsif($name eq "Kmer Content"){
				$self->kmer($pass);
			}elsif($name eq "Basic Statistics"){
				my $stats = sampStats->new();
				$stats->fillContainers(\@linesegs, $pass);
				$self->sampstats($stats);
			}

		}elsif($inloop){
			$store .= $line;
		}
	}

	close IN;
}

sub getOutArray{
	my ($self) = @_;
	#my @seHeaders = ("Sample", "Library", "FastqFile", "ReadNum", "TotalReads", "FilteredReads", "ReadLength", "TotalGCPerc", "BPQualityPass",
	#	"ReadQSPass", "QualityScorePass", "BPGCContentPass", "BPNContentPass", "OverRepPass", "KmerPass", "SeqDupPerc",
	#	"First25Qual", "Second25Qual", "Third25Qual", "Fourth25Qual", "#ReadsQS2-11", "#ReadsQS12-21", 
	#	"#ReadsQS22-31", "#ReadsQS32-40", "First25GC", "Second25GC", "Third25GC", "Fourth25GC", 
	#	"First25N", "Second25N", "Third25N", "Fourth25N" );
	
	my @output;
	
	push(@output, ($self->sample, $self->library, $self->file, $self->readNum, $self->sampstats->totseq, $self->sampstats->filtseq, $self->sampstats->seqlen, $self->sampstats->gc, $self->sampstats->pass));
	push(@output, ($self->pbpquality->pass, $self->pseqquality->pass, $self->pbpgccontent->pass, $self->pbpncontent->pass, $self->overrep, $self->kmer, $self->seqdup));
	my @keys = sort{ $a cmp $b} $self->pbpquality->getKeys();
	my @qskeys = sort{ $a cmp $b} $self->pseqquality->getKeys();
	push(@output, ($self->pbpquality->get($keys[0]), $self->pbpquality->get($keys[1]), $self->pbpquality->get($keys[2]), $self->pbpquality->get($keys[3])));
	push(@output, ($self->pseqquality->get($qskeys[0]), $self->pseqquality->get($qskeys[1]), $self->pseqquality->get($qskeys[2]), $self->pseqquality->get($qskeys[3])));
	push(@output, ($self->pbpgccontent->get($keys[0]), $self->pbpgccontent->get($keys[1]), $self->pbpgccontent->get($keys[2]), $self->pbpgccontent->get($keys[3])));
	push(@output, ($self->pbpncontent->get($keys[0]), $self->pbpncontent->get($keys[1]), $self->pbpncontent->get($keys[2]), $self->pbpncontent->get($keys[3])));
	
	return @output;
}

sub cleanUp{
	my ($self, $czip) = @_;
	
	if( -s $self->folder){
		$self->log->INFO("cleanup", "Removing fastqc uncompressed folder");
	}
	
	if( -s $self->zip && $czip){
		$self->log->INFO("cleanup", "Removing fastqc zip");
	}
}

sub getHeaderArray{
	my ($self) = @_;
	
	return @seHeaders;
}

__PACKAGE__->meta->make_immutable;

### All of the following packages are sub-data holders that process a large string of lines from fastqc text summary files

package sampStats;
use Mouse;
use namespace::autoclean;

has 'pass' => (is => 'rw', isa => 'Str', default => 'PASS');
has ['totseq', 'filtseq', 'gc'] => (is => 'rw', isa => 'Num', default => 0);
has 'seqlen' => (is => 'rw', isa => 'Any');

sub fillContainers{
	my ($self, $input, $pass) = @_;
	
	my @linesegs = @{$input};
	$self->pass($pass);
	for(my $x = 1; $x < scalar(@linesegs); $x++){
		if($linesegs[$x] =~ /^#/){next;}
		my @msegs = split(/\t/, $linesegs[$x]);
		if($msegs[0] eq "Total Sequences"){
			$self->totseq($msegs[1]);
		}elsif($msegs[0] eq "Filtered Sequences"){
			$self->filtseq($msegs[1]);
		}elsif($msegs[0] eq "Sequence length"){
			$self->seqlen($msegs[1]);
		}elsif($msegs[0] eq "%GC"){
			$self->gc($msegs[1]);
		}
	}
}
	

__PACKAGE__->meta->make_immutable;

package pbpQuality;
use Mouse;
use namespace::autoclean;
	
has 'pass' => (is => 'rw', isa => 'Str', default => 'PASS');
has 'quarts' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'getKeys' => 'keys',
		'exists' => 'exists',
		'get' => 'get',
		'set' => 'set',
	});

sub fillContainers{
	my ($self, $input, $pass) = @_;
	my @linesegs = @{$input};
	$self->pass($pass);
	
	my %meanhash;
	my %medianhash;
	for(my $x = 1; $x < scalar(@linesegs); $x++){
		if($linesegs[$x] =~ /^#/){next;}
		my @msegs = split(/\s+/, $linesegs[$x]);
		if($msegs[0] =~ /(\d+)-(\d+)/){
			for(my $y = $1; $y <= $2; $y++){
				$meanhash{$y} = $msegs[1];
				$medianhash{$y} = $msegs[2];
			}
		}else{
			$meanhash{$msegs[0]} = $msegs[1];
			$medianhash{$msegs[0]} = $msegs[2];
		}
	}
	my @values = sort {$b <=> $a} keys %meanhash;
	my %nmean;
	my %nmedian;
	for(my $x = 1; $x < $values[0] ; $x += int($values[0] * .25 + 0.5)){
		my @meantemp;
		my @mediantemp;
		for(my $y = $x; $y <= $values[0]; $y++){
			push(@meantemp, $meanhash{$y});
			push(@mediantemp, $medianhash{$y});
		}
		my $str = "$x-" . ($x + int($values[0] * .25 + 0.5));
		$self->set($str => $self->average(\@meantemp));
		#$nmedian{$str} = average(\@mediantemp);
	}

}

sub average{
	my $self = shift(@_);
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	if(scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}

__PACKAGE__->meta->make_immutable;

package pseqQuality;
use Mouse;
use namespace::autoclean;

has 'pass' => (is => 'rw', isa => 'Str', default => 'PASS');
has 'quarts' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'getKeys' => 'keys',
		'exists' => 'exists',
		'get' => 'get',
		'set' => 'set',
	});
	
sub fillContainers{
	my ($self, $input, $pass) = @_;
	my @linesegs = @{$input};
	$self->pass($pass);
	
	my %sqc;
	for(my $x = 1; $x < scalar(@linesegs); $x++){
		if($linesegs[$x] =~ /^#/){next;}
		my @msegs = split(/\s+/, $linesegs[$x]);
		$sqc{$msegs[0]} = $msegs[1];
	}
	my @values = sort {$b <=> $a} keys %sqc;
	my %final;
	for(my $x = 2; $x <= 40 +2; $x+=10){
		my @meanvalue;
		for(my $y = $x; $y <= $values[0]; $y++){
			push(@meanvalue, $sqc{$y}) if exists($sqc{$y});
		}
		my $end = ($x + 10 > 40)? 40 : $x + 10;
		my $str = "$x-$end";
		$self->set($str => $self->average(\@meanvalue));
		if($end == 40){last;}
	}
}

sub average{
	my $self = shift(@_);
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	if(scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}

__PACKAGE__->meta->make_immutable;

package pbpGCContent;
use Mouse;
use namespace::autoclean;

has 'pass' => (is => 'rw', isa => 'Str', default => 'PASS');
has 'quarts' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'getKeys' => 'keys',
		'exists' => 'exists',
		'get' => 'get',
		'set' => 'set',
	});

sub fillContainers{
	my ($self, $input, $pass) = @_;
	my @linesegs = @{$input};
	$self->pass($pass);
	
	my %meanhash;
	for(my $x = 1; $x < scalar(@linesegs); $x++){
		if($linesegs[$x] =~ /^#/){next;}
		my @msegs = split(/\s+/, $linesegs[$x]);
		if($msegs[0] =~ /(\d+)-(\d+)/){
			for(my $y = $1; $y <= $2; $y++){
				$meanhash{$y} = $msegs[1];
			}
		}else{
			$meanhash{$msegs[0]} = $msegs[1];
		}
	}
	my @values = sort {$b <=> $a} keys %meanhash;
	my $z = 1;
	for(my $x = 1; $x < $values[0]; $x += int($values[0] * .25 + 0.5)){
		my @meantemp;
		for(my $y = $x; $y <= $values[0]; $y++){
			push(@meantemp, $meanhash{$y});
		}
		my $str = "$x-" . ($x + int($values[0] * .25 + 0.5));
		$z++;
		$self->set($str => $self->average(\@meantemp));
	}
}

sub average{
	my $self = shift(@_);
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	if(scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}

__PACKAGE__->meta->make_immutable;

package pbpNContent;
use Mouse;
use namespace::autoclean;

has 'pass' => (is => 'rw', isa => 'Str', default => 'PASS');
has 'quarts' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'getKeys' => 'keys',
		'exists' => 'exists',
		'get' => 'get',
		'set' => 'set',
	});

sub fillContainers{
	my ($self, $input, $pass) = @_;
	my @linesegs = @{$input};
	$self->pass($pass);
	
	my %meanhash;
	for(my $x = 1; $x < scalar(@linesegs); $x++){
		if($linesegs[$x] =~ /^#/){next;}
		my @msegs = split(/\s+/, $linesegs[$x]);
		if($msegs[0] =~ /(\d+)-(\d+)/){
			for(my $y = $1; $y <= $2; $y++){
				$meanhash{$y} = $msegs[1];
			}
		}else{
			$meanhash{$msegs[0]} = $msegs[1];
		}
	}
	my @values = sort {$b <=> $a} keys %meanhash;
	my $z = 1;
	for(my $x = 1; $x < $values[0]; $x += int($values[0] * .25 + 0.5)){
		my @meantemp;
		for(my $y = $x; $y <= $values[0]; $y++){
			push(@meantemp, $meanhash{$y});
		}
		my $str = "$x-" . ($x + int($values[0] * .25 + 0.5));
		$z++;
		$self->set($str => $self->average(\@meantemp));
	}
}

sub average{
	my $self = shift(@_);
	my $array_r = shift(@_);
	my @numb = @{$array_r};
	if(scalar(@numb) == 0){
		return 0;
	}
	my $total3 = 0;
	foreach my $num1 (@numb) {
	$total3 += $num1;
	}
	my $mean3 = $total3 / (scalar @numb);
	return $mean3;
}

__PACKAGE__->meta->make_immutable;
1;
