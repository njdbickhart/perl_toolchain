#!/usr/bin/perl
# A list of utilities designed to process SAM/BAM files in my perl scripts
# Assumes SAMTOOLS is in the path!

package SamFileReader;
use Mouse;
use strict;
use File::Basename;
use namespace::autoclean;

has 'inputFile' => (is => 'rw', isa => 'SamFile');

# In case the program needs to create and index a Bam, the class has a container, here
has 'alternate' => (is => 'rw', isa => 'SamFile', predicate => 'has_bam');

# Samtools version checking wrapper
has 'samExe' => (is => 'rw', isa => 'SamtoolsExecutable', lazy => 1, default => sub{SamtoolsExecutable->new()});

# Initiator method to check file and determine if additional processing is needed
sub prepSam{
	my ($self, $file) = @_;
	$self->_checkFile($file);
}

# Uses samtools idxstats to get coverage data for a sam/bam
# RETURN:
# raw coverage, mapped coverage, hash->{chr} = [raw chr coverage, mapped chr coverage]
sub getXCov{
	my ($self) = @_;	
	
	my $SamFile = ($self->has_bam)? $self->alternate : $self->inputFile;

	my $readlen = $self->samExe->GetBamReadLen($SamFile->File);
	if(!$SamFile->isIndexed){
		$SamFile->checkIndex();
	}
	
	open(IN, $self->samExe->SamIdxstats() . $SamFile->File . " |") || die "[SAMFILEREADER] Could not generate index stats for bam file: " . $SamFile->File . "!\n";
	my $genomeLength = 0;
	my $mappedReadCount = 0;
	my $unmappedCount = 0;
	
	my %chrcov;
	
	while(my $line = <IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		my $mapchrcov = ($segs[2] * $readlen) / $segs[1];
		my $rawchrcov = (($segs[2] + $segs[3]) * $readlen) / $segs[1];
		$genomeLength += $segs[1];
		$mappedReadCount += $segs[2];
		$unmappedCount += $segs[3];
		$chrcov{$segs[0]} = [$rawchrcov, $mapchrcov];
	}
	close IN;
	
	my $rawCov = (($mappedReadCount + $unmappedCount) * $readlen)/ $genomeLength;
	my $mappedCov = ($mappedReadCount * $readlen) / $genomeLength;
	return ($rawCov, $mappedCov, \%chrcov);
}

sub _checkFile{
	my ($self, $file) = @_;
	my $sam = SamFile->new('File' => $file);
	$self->inputFile($sam);

	my $basename = basename($file);
	# Check if the file is a bam. If not, generate a new file that is a bam
	if($basename =~ m/\.bam$/ || -B $file){
		$sam->isBam(1);
		$sam->checkIndex;
	}else{
		# Generate a new bam file for use later
		my $bam = $sam->createBam();
		$self->alternate($bam);
	}
}


__PACKAGE__->meta->make_immutable;

package SamFile;
use Mouse;
use strict;
use File::Basename;
use namespace::autoclean;

has 'File' => (is => 'ro', isa => 'Str');
has 'isBam' => (is => 'rw', isa => 'Bool', default => 0);
has 'isIndexed' => (is => 'rw', isa => 'Bool', default => 0);
has 'samExe' => (is => 'rw', isa => 'SamtoolsExecutable', lazy => 1, default => sub{SamtoolsExecutable->new()});



# Wrapper to check if the bam is indexed
sub checkIndex{
	my ($self) = @_;
	
	if(-s $self->File . ".bai"){
		# Nothing to do, returning
		$self->isIndexed(1);
		#return;
	}else{
		print STDERR "[SAMFILE] Could not find index for bam, attempting to create one...\n";
		system($self->samExe->SamIndex() . $self->File);
		unless(-s $self->File . ".bai"){
			print STDERR "[SAMFILE] Could not create index! Perhaps the bam is not sorted?\n";
			print STDERR "[SAMFILE] premature exit...\n";
			exit;
		}
	}
}

# Returns a new SamFile object that is a BAM and is indexed
sub createBam{
	my ($self) = @_;
	my ($filename, $dirs, $suffix) = fileparse($self->File);
	
	my $bam = "$dirs/$filename.bam";
	print STDERR "[SAMFILE] Converting filetype: SAM -> " . $self->File . " to bam...\n";
	# Going to work with a default 4 threads for conversion
	system($self->samExe->SamToBam($self->File, $bam, 4));
	
	if(-s $bam){
		print STDERR "[SAMFILE] Successfully created bam file: $bam\n";
	}else{
		print STDERR "[SAMFILE] ERROR! Could not create bam file: $bam\n";
		print STDERR "[SAMFILE] premature exit...\n";
		exit;
	}
	
	print STDERR "[SAMFILE] Indexing bam file...\n";
	system($self->samExe->SamIndex() . $bam);
	
	if(-s "$bam.bai"){
		print STDERR "[SAMFILE] Successfully created index\n";
	}else{
		print STDERR "[SAMFILE] Could not create index: $bam.bai!\n";
		print STDERR "[SAMFILE] premature exit...\n";
		exit;
	}
	
	# Now to create and return the new SamFile object
	my $sam = SamFile->new('File' => $bam);
	$sam->isBam(1);
	$sam->isIndexed(1);
	return $sam;
}


__PACKAGE__->meta->make_immutable;

# This class is needed to cope with new HTSLIB options (and lack of backwards compatibility!)
package SamtoolsExecutable;
use Mouse;
#use StaticUtils;
use File::Basename;
use strict;
use namespace::autoclean;

has 'isHTSLib' => (is => 'rw', isa => 'Bool', lazy => 1, builder => '_checkVersion');

# Samples bam and gets read lengths from first reads
# Very simple implementation -- perhaps I could sample until all read groups are accounted for in the future?
sub GetBamReadLen{
	my ($self, $bam) = @_;
	
	open(IN, "samtools view $bam | head |") || die "[SAMEXE] Could not open bam for readlength checking!\n";
	my $line = <IN>;
	my @segs = split(/\t/, $line);
	close IN;
	
	return length($segs[9]);
}

# Returns a string to use for samtools idxstats
# Again, redundant, but sticks with the theme of the wrapper
sub SamIdxstats{
	my ($self) = @_;
	return "samtools idxstats ";
}

# Returns a program argument string that will convert a sam file to a bam file with sorting
# HTSlib allows threading of samtools sorting
sub SamToBam{
	my ($self, $sam, $bam, $threads) = @_;
	
	if($self->isHTSlib){
		return "samtools view -bS $sam | samtools sort -o $bam -T $bam.pre -\@ $threads - ";
	}else{
		my ($filename, $dirs, $suffix) = fileparse($bam);
		return "samtools view -bS $sam | samtools sort - $dirs/$filename ";
	}
}

# Since the syntax is the same, this is just redundant, but it sticks with the API format and is
# likely to not impede performance too much
sub SamIndex{
	my ($self) = @_;
	return "samtools index ";
}

sub _checkVersion{
	my ($self) = @_;
	checkReqs("samtools");
	
	# Now to grep out the version number:
	open(IN, "samtools | grep 'Version' |") || die "[SAMEXE] Could not open samtools for version checking!\n";
	while(my $line = <IN>){
		my ($num) = $line =~ /.*Version: (\d{1})\..+/;
		if($num == 1 || $num == 0){
			$self->isHTSLib($num);
		}
	}
	close IN;
}



__PACKAGE__->meta->make_immutable;

package StaticUtils;
use strict;

sub checkReqs{
	my ($program) = @_;
	
	my $found = 0;
	foreach my $path (split(/:/, $ENV{PATH})) {
		if( -f "$path/$program") {
			$found = 1;
			last;
		}
	}
	
	if(!$found){
		print "Error! Could not find the program: $program on your path! Please install it before continuing!\n";
		exit;
	}
}

# Returning a "one" for Perl module loading
1;
