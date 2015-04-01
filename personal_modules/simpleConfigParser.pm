#!/usr/bin/perl
# This module is designed to read a text file full of values delimited by equals signs
# Example: java=/my/path/to/java

package simpleConfigParser;
use strict;
use Mouse;
use namespace::autoclean;

has 'Map' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub {{}},
	handles => {
		'put' => 'set',
		'getKey' => 'get',
		'keyExists' => 'exists',});
		
sub loadConfigFile{
	my ($self, $file) = @_;
	open(IN, "< $file") || die "[CNFGPARSE] Could not open config file: $file!\n";
	while(my $line = <IN>){
		if($line =~ /^#/){
			next;
		}
		chomp $line;
		my @segs = split(/\=/, $line);
		$self->put($segs[0] => $segs[1]);
	}
	close IN;
}

sub checkReqKeys{
	my ($self, $keys) = @_;
	my @missing;
	foreach my $k (@{$keys}){
		if(!$self->keyExists($k)){
			push(@missing, $k);
		}
	}
	if(scalar(@missing) > 0){
		print STDERR "[CNFGPARSE] Error! Missing required config values:\n";
		foreach my $m (@missing){
			print STDERR "\t$m\n";
		}
		exit;
	}
}

__PACKAGE__->meta->make_immutable;

1;