Instructions for Compiling PEST++ on Linux
==========================================

Introduction
------------

The Windows version of PEST++ supports the serial, YAMR and GENIE run managers.
The “Updated PEST++ input instructions” located at
[`http://www.pestpp.org/download/pest++_input.pdf`] provides a quick summary of
each as well as the command line arguments necessary to invoke them.  A Linux
version Genie is not available, the Linux version of PEST++ only supports the
serial and YAMR run managers.

Dependencies:

1. gcc, g++ and gfortran 4.9.0 or newer : PEST++ relies on some newer
features of C++11 standard which are only supported in these newer releases .
2. blas
3. lapack

Instructions to build PEST++ on Linux:

1. Edit the file `./pest++/common/config_os.h` and make the following
change:

change:

    //#define OS_LINUX
    #define OS_WIN

to:

    #define OS_LINUX
    //#define OS_WIN

2. In the main pestpp directory execute the following command to delete any files associated with a previous build:

    make –f makefile_linux clean

3. Build PEST++ in the main pestpp directory using the following command:

    make –f makefile_linux clean

This creates the PEST++ executable:

    pestpp/pest++/pestpp 

