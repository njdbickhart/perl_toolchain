#!/bin/bash
# This is the initial installation script that downloads a version of local::lib and uses bootstrapping to install perl modules

echo "Downloading local::lib for subsequent bootstrapping..."
echo "Assuming that you are running this script from inside the perl_toolchain folder in your HOME directory! If not, please use control-C to cancel this install!"

cd ../

echo 'wget http://search.cpan.org/CPAN/authors/id/H/HA/HAARG/local-lib-2.000015.tar.gz'
wget http://search.cpan.org/CPAN/authors/id/H/HA/HAARG/local-lib-2.000015.tar.gz

echo 'tar -xvf local-lib-2.000015.tar.gz'
tar -xvf local-lib-2.000015.tar.gz
mkdir tar_repo
mv local-lib-2.000015.tar.gz ./tar_repo/

cd local-lib-2.000015/

echo "Beginning local::lib bootstrapping..."
perl Makefile.PL --bootstrap
make test
make install
echo '[ $SHLVL -eq 1 ] && eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >>~/.bashrc
echo 'export PERL5LIB=$PERL5LIB:$HOME/perl_toolchain/personal_modules/' >>~/.bashrc


cd ../

if [ ! -d "$HOME/bin" ]; then
  mkdir $HOME/bin
fi

echo "Installing cpanm..."
wget http://xrl.us/cpanm
chmod +x ./cpanm
mv ./cpanm $HOME/bin/

echo "Finished! Please either restart your system or relog to apply changes."