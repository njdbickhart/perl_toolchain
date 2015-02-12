#!/usr/bin/perl
# This is a collection of handy utilities designed to make memory management and time tracking easier

package memTester;
use Mouse;
use strict;
use namespace::autoclean;

has 'maxMem' => (is => 'rw'

__PACKAGE__->meta->make_immutable;
