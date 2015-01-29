#!/bin/bash
# This script uses cpanm to queue the installation of key Perl modules into a local repository for later use

echo "Installing necessary Perl libraries..."

cpanm namespace::autoclean
cpanm Mouse
cpanm Sys::CpuLoadX
cpanm Sys::CpuAffinity

cpanm MouseX::NativeTraits
cpanm Spreadsheet::WriteExcel
cpanm List::Util
echo "Attempting Forks Super install..."
echo "WARNING! This module is stubborn and may need to be installed twice."
echo "To do this, run '$ cpanm Forks::Super' after the termination of this script."
cpanm --verbose Forks::Super


echo "Finished installation!"
