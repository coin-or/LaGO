# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([ThirdPartyRANLIB],[0.0],[http://projects.coin-or.org/LaGO/newticket])

AC_COPYRIGHT([
Copyright 2007 Stefan Vigerske
All Rights Reserved.
This file is part of the open source package LaGO which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(ranlib.c/src/ranlib.h)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(FILIB)

# We need egrep and ln -s somewhere
#AC_PROG_EGREP
#AC_PROG_LN_S

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Get the name of the C++ compiler and appropriate compiler options
#AC_COIN_PROG_CXX

# Get the C++ runtime libraries in case we want to link with a C or
# Fortran compiler
# AC_COIN_CXXLIBS

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS


############
# END PART #
############

AC_CONFIG_FILES([Makefile])

AC_COIN_FINALIZE
