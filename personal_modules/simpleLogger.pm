#!/usr/bin/perl
# This is a module designed to provide simple logging functionality to Perl

package simpleLogger;
use strict;
use Mouse;
use namespace::autoclean;
use IO::File;
use Fcntl qw(:flock SEEK_END);

has 'logFileBaseStr' => (is => 'ro', isa => 'Str', required => 1);
has 'logHandle' => (is => 'rw', isa => 'IO::File');

sub OpenLogger{
	my ($self, $outputfolder) = @_;
	my $epoc = time();
	my $file = "$outputfolder/" . $self->logFileBaseStr() . ".$epoc.log";
	#open(LOG, "> $file") || die "[SMPLOG] Could not open log file!\n";
	$self->logHandle(IO::File->new());
	$self->logHandle->open("> $file");
}

sub GetDate{
	my @datesegs = localtime();
	my $dateformat = sprintf("%02d/%02d/%02d - %02d:%02d:%02d", $datesegs[3], $datesegs[4], $datesegs[5], $datesegs[2], $datesegs[1], $datesegs[0]);
	return $dateformat;
}

sub Info{
	my ($self, $position, $message) = @_;
	my $time = $self->GetDate;
	flock($self->logHandle, LOCK_EX);
	print {$self->logHandle} "$time \| $position - $message\n";
	$self->logHandle->flush();
	flock($self->logHandle, LOCK_UN);
}

sub Fatal{
	my ($self, $position, $message) = @_;
	my $time = $self->GetDate;
	flock($self->logHandle, LOCK_EX);
	print STDERR "Error with pipeline! $position - $message\n";
	print {$self->logHandle} "FATAL! $time \| $position - $message\n";
	$self->logHandle->flush();
	flock($self->logHandle, LOCK_UN);
	exit;
}

sub Warn{
	my ($self, $position, $message) = @_;
	my $time = $self->GetDate;
	flock($self->logHandle, LOCK_EX);
	print STDERR "Pipeline warning! $position - $message\n";
	print {$self->logHandle} "WARNING! $time \| $position - $message\n";
	$self->logHandle->flush();
	flock($self->logHandle, LOCK_UN);
}

sub Close{
	my ($self) = @_;
	$self->logHandle->close();
}

__PACKAGE__->meta->make_immutable;

1;