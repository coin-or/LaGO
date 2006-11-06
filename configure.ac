# Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske
# Template from: Andreas Waechter

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([LaGO],[0.2],[http://projects.coin-or.org/LaGO/newticket])

AC_COPYRIGHT([
Copyright 2006 Ivo Nowak and Stefan Vigerske
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

#AC_REVISION(0.9)

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(configure.ac)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#     We only need automake to generate Makefiles for the distribution      #
#############################################################################

# Initialize automake
AC_COIN_INIT_AUTOMAKE

#############################################################################
#             Check which ThirdParty sources are available                  #
#############################################################################

AC_COIN_THIRDPARTY_SUBDIRS([ASL],[ThirdParty/ASL],[solvers/asl.h])
#AC_COIN_THIRDPARTY_SUBDIRS([TNT],[ThirdParty/TNT],[tnt.h])
#AC_COIN_THIRDPARTY_SUBDIRS([RANLIB],[ThirdParty/ranlib],[ranlib.c/src/ranlib.h])
#AC_COIN_THIRDPARTY_SUBDIRS([METIS],[ThirdParty/METIS],[metis-4.0/Lib/metis.h])
#AC_COIN_THIRDPARTY_SUBDIRS([GAMSIO],[ThirdParty/GAMS_IO],[GAMS_IO/iolib.h])

#############################################################################
#                  Check which subprojects are there                        #
#############################################################################

AC_COIN_MAIN_SUBDIRS(CoinUtils)
AC_COIN_MAIN_SUBDIRS(Clp)
AC_COIN_MAIN_SUBDIRS(Osi)
AC_COIN_MAIN_SUBDIRS(Cgl)
AC_COIN_MAIN_SUBDIRS(Ipopt)
AC_COIN_MAIN_SUBDIRS(LaGO)


##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE