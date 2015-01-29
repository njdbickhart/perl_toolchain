# perl_toolchain
Backup for my personal perl scripts

In order to install prerequisite modules to run these scripts, please use the installation shell scripts in the base directory of the repository in the order that they are listed.

The first step is to generate a local::lib bootstrap that will create the necessary folders and environmental variables for Perl to use for installing future modules. In order to accomplish this, you simply need to run the following command within the "perl_toolchain" directory in your $HOME folder:

```bash
$ sh step1_locallib_install.sh
```

After this completes, you need to restart your shell session (or reboot the machine) in order to apply the changes to your current session. After rebooting, you will then download the necessary Perl modules by using this command within the "perl_toolchain" folder:

```bash
$ sh step2_required_library_install.sh
```

This will take some time, but after it is finished, you should be able to run just about all of the scripts


# Troubleshooting
There is an issue with the Forks::Super library that might cause errors during the installation. If you receive an error at the end of the "step 2" script, please try reinstalling the module by running the following command:

```bash
$ cpanm --verbose Forks::Super
```

This should install the library properly after a lengthy series of tests.
