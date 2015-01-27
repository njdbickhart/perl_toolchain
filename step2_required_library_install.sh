#!/bin/bash
# This script uses cpanm to queue the installation of key Perl modules into a local repository for later use

echo "Installing necessary Perl libraries..."

cpanm namespace::autoclean
cpanm Mouse
cpanm Sys::CpuLoadX
cpanm Sys::CpuAffinity
cpanm --verbose Forks::Super
cpanm MouseX::NativeTraits
cpanm Spreadsheet::WriteExcel
cpanm List::Util

echo "Finished installation!"