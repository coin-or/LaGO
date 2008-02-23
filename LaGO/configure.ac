# Copyright (C) 2007-2008 Stefan Vigerske
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

AC_INIT([LaGO],[0.4],[http://projects.coin-or.org/LaGO/newticket])

AC_COPYRIGHT([
Copyright 2007-2008 Stefan Vigerske
All Rights Reserved.
This file is part of the open source package LaGO which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/Base/LaGOConfig.h)

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
AC_COIN_DEBUG_COMPILE(LaGO)

# We need egrep and ln -s somewhere
#AC_PROG_EGREP
#AC_PROG_LN_S

# disable warning messages to make LaGO compiling
#coin_skip_warn_cxxflags=yes

# Get the name of the C compiler and appropriate compiler options
AC_COIN_PROG_CC

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Get the C++ runtime libraries in case we want to link with a C or
# Fortran compiler
# AC_COIN_CXXLIBS

# Get the name of the Fortran compiler and appropriate compiler options
#AC_COIN_PROG_F77

# Find out how to call Fortran from C and determine Fortran runtime libraries
#AC_COIN_F77_WRAPPERS

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(LaGO)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(Osi)
# later: this one is optional, if SNOPT is there
AC_COIN_HAS_PROJECT(Ipopt)
# later: this one is optional, if, e.g., CPLEX is there
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(Cbc)
AC_COIN_HAS_PROJECT(Bonmin)

############################################################################
############################################################################
#                   Stuff that we need for C++ programs                    #
############################################################################
############################################################################

AC_LANG_PUSH(C++)

################
# Header Files #
################


AC_COIN_CHECK_CXX_CHEADER(math)
AC_COIN_CHECK_CXX_CHEADER(ctype)
AC_COIN_CHECK_CXX_CHEADER(assert)
AC_COIN_CHECK_CXX_CHEADER(stdio)
AC_COIN_CHECK_CXX_CHEADER(stdlib)
AC_COIN_CHECK_CXX_CHEADER(stdarg)
AC_COIN_CHECK_CXX_CHEADER(float)
AC_COIN_CHECK_CXX_CHEADER(time)

###################
# Function finite #
###################

# First check for a few headers

case $build in
  *-cygwin*) ;;
  *) AC_COIN_CHECK_CXX_CHEADER([ieeefp]) ;;
esac

# now try to guess the correct name of the finite function

MY_C_FINITE=
AC_CHECK_DECL([finite],[MY_C_FINITE=finite],,[
#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# endif
#endif
#ifdef HAVE_CFLOAT
# include <cfloat>
#else
# ifdef HAVE_FLOAT_H
#  include <float.h>
# endif
#endif
#ifdef HAVE_CIEEEFP
# include <cieeefp>
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
# endif
#endif])
if test -z "$MY_C_FINITE"; then
  AC_CHECK_DECL([_finite],[MY_C_FINITE=_finite],,[
#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# endif
#endif
#ifdef HAVE_CFLOAT
# include <cfloat>
#else
# ifdef HAVE_FLOAT_H
#  include <float.h>
# endif
#endif
#ifdef HAVE_CIEEEFP
# include <cieeefp>
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
# endif
#endif])
  if test -z "$MY_C_FINITE"; then
    AC_CHECK_DECL([isfinite],[MY_C_FINITE=isfinite],,[
#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# endif
#endif
#ifdef HAVE_CFLOAT
# include <cfloat>
#else
# ifdef HAVE_FLOAT_H
#  include <float.h>
# endif
#endif
#ifdef HAVE_CIEEEFP
# include <cieeefp>
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
# endif
#endif])
  fi
fi
if test -z "$MY_C_FINITE"; then
  AC_MSG_WARN(Cannot find C-function for checking NaN and Inf.  Check will be disabled)
else
  AC_DEFINE_UNQUOTED([MY_C_FINITE],[$MY_C_FINITE],
                     [Define to be the name of C-function for NaNInf check])
fi

##########################################################################

# Additional libraries and packages

# We need to restore the value of LIBS before those test, otherwise libtool
# will hick-up
LIBS_without_ADDLIBS="$LIBS"

##################
# LAPACK Library #
##################

#AC_COIN_HAS_LAPACK

#################
# CPLEX Library #
#################

AC_LANG_PUSH(C)
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_LANG_POP(C)

####################
# RANLIB.C LIBRARY #
####################

coin_ranlibdir=../ThirdParty/RANLIB.C
coin_ranlibsrcdir=`$CYGPATH_W $abs_source_dir/$coin_ranlibdir | sed -e sX\\\\\\\\X/Xg`

AC_ARG_WITH([ranlib],
            AC_HELP_STRING([--with-ranlib-dir],
                           [specify directory of Random Number generator library libranlib.c.a]),
            [ranlib_dir=$withval],
	          [ranlib_dir=BUILD])

AC_LANG_PUSH(C)						
						
if test x"$ranlib_dir" = x"BUILD"; then
	AC_CHECK_FILE([$coin_ranlibsrcdir/ranlib.c/src/ranlib.h],,
		AC_MSG_ERROR([sources for library RANLIB.C not available at $coin_ranlibsrcdir]))
	RANLIBC_LIB="../../$coin_ranlibdir/libranlib.c.la -lm"
	RANLIBC_CPPFLAGS="-I$coin_ranlibsrcdir/ranlib.c/src"
else
  AC_MSG_CHECKING([whether ranlib.c library at $ranlib_dir works])

	RANLIBC_LIB="$ranlib_dir/libranlib.c.a -lm"
	RANLIBC_CPPFLAGS="-I$ranlib_dir/src"
	LIBS_save="$LIBS"
	LIBS="$RANLIBC_LIB $LIBS"
	CPPFLAGS="$RANLIBC_CPPFLAGS $CPPFLAGS"
	AC_TRY_LINK([#include "ranlib.h"], [genunf(0.,1.)],
							[AC_MSG_RESULT([yes])],
							[AC_MSG_RESULT([no])
							AC_MSG_ERROR([ranlib.c at $use_ranlib does not work])])
	LIBS="$LIBS_save"
fi

AC_SUBST(RANLIBC_LIB)
AC_SUBST(RANLIBC_CPPFLAGS)

AC_LANG_POP(C)


###########
# FILIB++ #
###########

#use path relative to directory LaGO/src/Applications where linking to the filib is done
coin_filibdir=../ThirdParty/FILIB++
coin_filibsrcdir=`$CYGPATH_W $abs_source_dir/$coin_filibdir | sed -e sX\\\\\\\\X/Xg`

AC_ARG_WITH([filib],
            AC_HELP_STRING([--with-filib],
                           [specify directory of FILIB++]),
            [filib_dir=$withval],
            [filib_dir=BUILD])

if test x"$filib_dir" = x"BUILD"; then
	AC_CHECK_FILE([$coin_filibsrcdir/interval/interval/interval.hpp],
    [coin_has_filib="yes"
     FILIB_LIB="../../$coin_filibdir/interval/libprim.la"
     FILIB_CPPFLAGS="-I$coin_filibsrcdir/interval"],
    [coin_has_filib="no"])
else
  AC_CHECK_FILE([$filib_dir/interval/interval/interval.hpp],
    [coin_has_filib="yes"
     FILIB_LIB="$filib_dir/lib/libprim.la"
     FILIB_CPPFLAGS="-I$filib_dir/interval"],
    [coin_has_filib="no"])
#TODO: check also whether library is working
fi

if test x"$coin_has_filib" = xyes; then
	AC_DEFINE([COIN_HAS_FILIB],[1],[If defined, the FILIB++ is available.])
	AC_SUBST(FILIB_LIB)
	AC_SUBST(FILIB_CPPFLAGS)
fi
AM_CONDITIONAL([COIN_HAS_FILIB],[test "$coin_has_filib" = yes])

####################################
# GAMS SYSTEM and GAMS I/O LIBRARY #
####################################

AC_COIN_HAVE_GAMS

#if we have gams i/o libraries, then use them
if test $coin_has_gamsio = yes ; then
  AC_COIN_USE_GAMS
  if test x$coin_has_cpx = xtrue; then
    AC_COIN_GAMSCPLEX_LICENCE
  else
    AM_CONDITIONAL([COIN_HAS_GAMSCPLEXLICENCE],[test 0 = 1])
  fi
else
  AM_CONDITIONAL([COIN_HAS_GAMSIO],[test 0 = 1])
  AM_CONDITIONAL([COIN_HAS_GAMSCPLEXLICENCE],[test 0 = 1])
fi

#######################
# AMPL SOLVER LIBRARY #
#######################

#LIBS="$LIBS_without_ADDLIBS"

#if test "$coin_has_gamsio" = "no"; then
#  AC_COIN_HAS_ASL
	
#	if test "$coin_has_asl" = "no"; then
#		AC_MSG_ERROR([Either GAMS I/O libraries or AMPL Solver Library (ASL) needed for LaGO])
#	fi
	
#	CPPFLAGS="$ASL_CPPFLAGS $CPPFLAGS"
#	ADDLIBS="$ASLLIB $ADDLIBS"
#else
#  AC_MSG_NOTICE([skip test for AMPL solver library because GAMS I/O libraries are working])

#	AM_CONDITIONAL([COIN_HAS_ASL],[test a=b])
#fi

############
# END PART #
############

AC_LANG_POP(C++)

AC_COIN_VPATH_LINK(src/Applications/gmsla_us.run)
AC_COIN_VPATH_LINK(src/Applications/install/gmsla_un.dat)

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/Base/Makefile
                 src/Problem/Makefile
                 src/Algorithms/Makefile
                 src/Interfaces/Makefile
                 src/Interfaces/GAMS/Makefile
                 src/Applications/Makefile
                ])

AC_CONFIG_HEADER([inc/config_lago.h])

AC_COIN_FINALIZE
