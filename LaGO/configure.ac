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

AC_INIT([LaGO],[0.3],[http://projects.coin-or.org/LaGO/newticket])

AC_COPYRIGHT([
Copyright 2006 Ivo Nowak and Stefan Vigerske
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/standard.h)

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
coin_skip_warn_cxxflags=yes

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
AC_COIN_F77_WRAPPERS

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

# A Fix for the AIX compiler
#AM_CONDITIONAL([AIX_FIX_COMPILER_BUG], test "x$AIX_FIX_COMPILER_BUG" != x)

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


#################
# CPLEX Library #
#################

AC_LANG_PUSH(C)
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_LANG_POP(C)


###############
# TNT Library #
###############

coin_tntdir=../ThirdParty/TNT
coin_tntsrcdir=$abs_source_dir/$coin_tntdir

AC_ARG_WITH([tnt],
            AC_HELP_STRING([--with-tnt],
                           [specify directory of Template Numerical Toolkit (TNT) header files]),
            [use_tnt=$withval],
	          [use_tnt=$coin_tntsrcdir])

AC_CHECK_FILE([$use_tnt/tnt.h],,AC_MSG_ERROR([Template Numerical Toolkit (TNT) header files not available]))

CPPFLAGS="-I$use_tnt $CPPFLAGS"

####################
# RANLIB.C LIBRARY #
####################

coin_ranlibdir=../ThirdParty/ranlib/ranlib.c
coin_ranlibsrcdir=$abs_source_dir/$coin_ranlibdir

AC_ARG_WITH([ranlib],
            AC_HELP_STRING([--with-ranlib],
                           [specify directory of Random Number generator library ranlib.c]),
            [use_ranlib=$withval],
	          [use_ranlib=$coin_ranlibsrcdir])

AC_LANG_PUSH(C)						
						
AC_MSG_CHECKING([whether ranlib.c library at $use_ranlib works])
LIBS="-L$use_ranlib -lranlib -lm $LIBS"
ADDLIBS="-L$use_ranlib -lranlib -lm $ADDLIBS"
CPPFLAGS="-I$use_ranlib/src $CPPFLAGS"
AC_TRY_LINK([#include "ranlib.h"], [genunf(0.,1.)],
            [AC_MSG_RESULT([yes])],
            [AC_MSG_RESULT([no])
             AC_MSG_ERROR([ranlib.c at $use_ranlib does not work])])

AC_LANG_POP(C)

#################
# METIS Library #
#################

coin_metisdir=../ThirdParty/METIS/metis-4.0
coin_metissrcdir=$abs_source_dir/$coin_metisdir

AC_ARG_WITH([metis],
            AC_HELP_STRING([--with-metis],
                           [specify directory of METIS]),
            [use_metis=$withval],
	          [use_metis=$coin_metissrcdir])

AC_LANG_PUSH(C)						
						
AC_MSG_CHECKING([whether METIS library at $use_metis works])
LIBS="-L$use_metis -lmetis $LIBS"
ADDLIBS="-L$use_metis -lmetis -lm $ADDLIBS"
CPPFLAGS="-I$use_metis/Lib $CPPFLAGS"
AC_TRY_LINK([#include "metis.h"], [idxtype a],
            [AC_MSG_RESULT([yes])],
            [AC_MSG_RESULT([no])
             AC_MSG_ERROR([METIS at $use_ranlib does not work])])

AC_LANG_POP(C)

###########
# FILIB++ #
###########

coin_filibdir=../ThirdParty/FILIB++
coin_filibsrcdir=$abs_source_dir/$coin_filibdir

AC_ARG_WITH([filib],
            AC_HELP_STRING([--with-filib],
                           [specify directory of FILIB++]),
            [use_filib=$withval
             filib_libdir=$use_filib/lib],
	          [use_filib=$coin_filibsrcdir
             filib_libdir=../$coin_filibdir/interval/.libs
            ])
	          
	          
AC_CHECK_FILE([$use_filib/interval/interval/interval.hpp],
               [AC_DEFINE([COIN_HAS_FILIB],[1], [If defined, the FILIB++ is available.])
     	          coin_has_filib=yes
                ADDLIBS="-L$filib_libdir -lprim $ADDLIBS"
                CPPFLAGS="-I$use_filib/interval $CPPFLAGS"],
               [coin_has_filib=no])

#AC_MSG_CHECKING([whether FILIB++ library at $use_filib works])
#coin_save_libs=$LIBS
#coin_save_cppflags=$CPPFLAGS
#LIBS="-L$filib_libdir -lprim $LIBS"
#CPPFLAGS="-I$use_filib/interval $CPPFLAGS"
#AC_TRY_LINK([#include "interval/interval.hpp"], [filib::fp_traits<double>::setup();],
#  	        [AC_MSG_RESULT([yes])
#  	         AC_DEFINE([COIN_HAS_FILIB],[1], [If defined, the FILIB++ is available.])
#     	       coin_has_filib=yes
#             ADDLIBS="-L$filib_libdir -lprim $ADDLIBS"
#       	    ],
#         	  [AC_MSG_RESULT([no])
#           	 coin_has_filib=no
#           	 CPPFLAGS=$coin_save_cppflags
#           	])
#LIBS=$coin_save_libs

AM_CONDITIONAL([COIN_HAS_FILIB],[test "$coin_has_filib" = yes])

###############
# GAMS SYSTEM #
###############

AC_ARG_WITH([gamssystem],
  AC_HELP_STRING([--with-gamssystem],[specify directory of GAMS distribution]),
  AC_CHECK_FILE([$withval/gams],[use_gamssys=$withval],[use_gamssys=]),
  [AC_PATH_PROG(gamspath, [gams],,,)
   use_gamssys=${gamspath/%gams/}
  ])
GAMS_PATH=$use_gamssys

####################
# GAMS I/O LIBRARY #
####################

# This macro checks for a library containing the GAMS I/O library.
# It defines the makefile conditional and preprocessor macro
# COIN_HAS_GAMSIO, if GAMS I/O is available.

GAMSIO_CODE=unsupportedarchitecture

AC_ARG_WITH([gamsio-code],
  AC_HELP_STRING([--with-gamsio-code],[specify the GAMS I/O code used to identify the GAMS I/O libraries types]),
  [GAMSIO_CODE=$withval],
  [case $build in
      x86_64-*-linux-*)
        GAMSIO_CODE=LEI
      ;;
      i?86-*-linux-*)
        GAMSIO_CODE=LX3
      ;;
      *-cygwin* | *-mingw*)
        GAMSIO_CODE=VIS
      ;;
      i?86-*-darwin*)
        GAMSIO_CODE=DII
      ;;
      i?86-*-solaris*)
        GAMSIO_CODE=SIG
      ;;
      *)
        AC_MSG_WARN([Build type $build not supported by GAMS I/O libraries.])
  esac]
)

case $GAMSIO_CODE in
  LEI)
      gamsio_system=Linux
    ;;
  LX3)
      gamsio_system=Linux
    ;;
  VIS)
      gamsio_system=Windows
    ;;
  DII)
      gamsio_system=Darwin
    ;;
  SIG)
      gamsio_system=Solaris
    ;;
  *)
      AC_MSG_WARN([GAMS I/O libraries with system code $GAMSIO_CODE not supported.])
esac

#Intel Fortran libraries (required for Linux systems when not compiling with Intel compiler)
IFORT_LIBS=
if test x"$gamsio_system" = xLinux; then
	case "$CXX" in
		icpc* | */icpc*)
			;;
		*)
		  AC_ARG_WITH([ifort-libdir],
        AC_HELP_STRING([--with-ifort-libdir],[specify directory of Intel Fortran libraries]),
        AC_CHECK_FILE([$withval/libifcore.so],[IFORT_LIBS="-L$withval"],),
        [if test x"$use_gamssys" != x; then
          AC_CHECK_FILE([$use_gamssys/libifcore.so],[IFORT_LIBS="-L$use_gamssys"],)
         else
          AC_MSG_WARN([Path for Intel Fortran Libraries not given by --with-ifort-libdir and no GAMS system available. GAMS I/O libraries might fail to link.])
         fi
        ])
      IFORT_LIBS="$IFORT_LIBS -lifcore -limf -lirc"
	esac
fi

coin_gamsioobjdir=../ThirdParty/GAMSIO/$GAMSIO_CODE
coin_gamsiosrcdir=`$CYGPATH_W $abs_source_dir/$coin_gamsioobjdir | sed -e sX\\\\\\\\X/Xg`

#if user provided GAMS I/O libraries, use this one, otherwise check ThirdParty/GAMSIO directory
AC_ARG_WITH([gamsio],
  AC_HELP_STRING([--with-gamsio],[specify directory of GAMS I/O libraries and header files]),
  [use_gamsio=$withval],
  [AC_CHECK_FILE([$coin_gamsiosrcdir/iolib.h],[use_gamsio=$coin_gamsiosrcdir],[use_gamsio=])])

if test x"$use_gamsio" = x; then
  coin_has_gamsio=no
  coin_has_gamscplexlice=no
	AC_MSG_WARN([no GAMS I/O libraries found.])
else  # Check if GAMS I/O libs work
	AC_MSG_CHECKING([whether GAMS I/O libraries at $use_gamsio work])
	if test x"$GAMSIO_CODE" = xLX3; then
	  GAMSIO_LIBS="$use_gamsio/iolib.a $use_gamsio/nliolib.a $use_gamsio/clicelib.a $use_gamsio/gclib.a -ldl $IFORT_LIBS"
	  SMAG_LIBS="$use_gamsio/clicelib.a $use_gamsio/libsmag.a $use_gamsio/gclib.a $use_gamsio/libg2d.a -ldl $IFORT_LIBS"
	elif test x"$GAMSIO_CODE" = xLEI; then
	  GAMSIO_LIBS="$use_gamsio/iolib.a $use_gamsio/nliolib.a $use_gamsio/clicelib.a $use_gamsio/gclib.a -ldl $IFORT_LIBS"
	  SMAG_LIBS="$use_gamsio/clicelib.a $use_gamsio/libsmag.a $use_gamsio/gclib.a $use_gamsio/libg2d.a -ldl $IFORT_LIBS"
	elif test x"$gamsio_system" = xWindows; then
	  GAMSIO_LIBS="$use_gamsio/iolib.lib $use_gamsio/nliolib.lib $use_gamsio/clicelib.lib $use_gamsio/gclib.lib"
	  SMAG_LIBS="$use_gamsio/clicelib.lib $use_gamsio/smag.lib $use_gamsio/gclib.lib $use_gamsio/g2d.lib"
	elif test x"$GAMSIO_CODE" = xDII; then
	  GAMSIO_LIBS="$use_gamsio/iolib.a $use_gamsio/nliolib.a $use_gamsio/clicelib.a $use_gamsio/gclib.a -ldl -lSystemStubs"
	  SMAG_LIBS="$use_gamsio/clicelib.a $use_gamsio/libsmag.a $use_gamsio/gclib.a $use_gamsio/libg2d.a -ldl -lSystemStubs"
	elif test x"$GAMSIO_CODE" = xSIG; then
	  GAMSIO_LIBS="$use_gamsio/iolib.a $use_gamsio/nliolib.a $use_gamsio/clicelib.a $use_gamsio/gclib.a"
	  SMAG_LIBS="$use_gamsio/clicelib.a $use_gamsio/libsmag.a $use_gamsio/gclib.a $use_gamsio/libg2d.a"
	fi
	GAMSIO_CPPFLAGS="-I$use_gamsio"
	LIBS="$GAMSIO_LIBS $LIBS"
	CPPFLAGS_save=$CPPFLAGS
	CPPFLAGS="$GAMSIO_CPPFLAGS $CPPFLAGS"
	AC_TRY_LINK([
#include <cstdio>
#include "iolib.h"], [gfinit()],
  	[AC_MSG_RESULT([yes])
   	LIBS=$LIBS_without_ADDLIBS
   	CPPFLAGS=$CPPFLAGS_save],
  	[AC_MSG_RESULT([no])
   	AC_MSG_ERROR([GAMS I/O library $use_gamsio does not work])])


	coin_has_gamsio=yes
	AC_DEFINE([COIN_HAS_GAMSIO],[1],[If defined, the GAMS I/O Library is available.])
	AC_SUBST(coin_has_gamsio)
	AC_SUBST(GAMS_PATH)
	AC_SUBST(SMAG_LIBS)
	AC_SUBST(GAMSIO_CODE)

	#check for libgamscplexlice
	if test x"$coin_has_gamsio" = xyes; then
		if test x"$coin_has_cpx" = xtrue; then
			AC_LANG_PUSH(C)						
	
			AC_MSG_CHECKING([whether GAMS/CPLEX licence library is present and working])
			coin_save_libs=$LIBS
			coin_save_cppflags=$CPPFLAGS
			CPPFLAGS="$GAMSIO_CPPFLAGS $CPPFLAGS"
			LIBS="$use_gamsio/libgamscplexlice.a $LIBS $CPXLIB $GAMSIO_LIBS"
			AC_TRY_LINK([#include "gamscplexlice.h"], [gamscplexlice(0,0,0,0,0,0,0,0,0,0,0,0,0)],
	  		          [AC_MSG_RESULT([yes])
        	         coin_has_gamscplexlice=yes
        	         GAMSIO_LIBS="$use_gamsio/libgamscplexlice.a $GAMSIO_LIBS"
		        	     AC_DEFINE([COIN_HAS_GAMSCPLEXLICE],[1],[If defined, the GAMS/CPLEX licence library is available.])
              	  ],
      	        	[AC_MSG_RESULT([no])
                 	coin_has_gamscplexlice=no]
                 );
			LIBS=$coin_save_libs
			CPPFLAGS=$coin_save_cppflags
			AC_LANG_POP(C)
		else
			coin_has_gamscplexlice=no
		fi
	else
		coin_has_gamscplexlice=no
	fi

	AC_SUBST(GAMSIO_CPPFLAGS)
	AC_SUBST(GAMSIO_LIBS)
fi

AM_CONDITIONAL([COIN_HAS_GAMSIO],[test x"$use_gamsio" != x])
AM_CONDITIONAL([COIN_HAS_GAMSCPLEXLICE],[test "$coin_has_gamscplexlice" = yes])
AM_CONDITIONAL([GAMSIO_IS_UNIX],[test x"$gamsio_system" != xWindows])
AM_CONDITIONAL([GAMSIO_IS_WINDOWS],[test x"$gamsio_system" = xWindows])

#######################
# AMPL SOLVER LIBRARY #
#######################

LIBS="$LIBS_without_ADDLIBS"

if test "$coin_has_gamsio" = "no"; then
  AC_COIN_HAS_ASL
	
	if test "$coin_has_asl" = "no"; then
		AC_MSG_ERROR([Either GAMS I/O libraries or AMPL Solver Library (ASL) needed for LaGO])
	fi
	
	CPPFLAGS="$ASL_CPPFLAGS $CPPFLAGS"
#	LIBS="$ASLLIB $LIBS"
	ADDLIBS="$ASLLIB $ADDLIBS"
else
  AC_MSG_NOTICE([skip test for AMPL solver library because GAMS I/O libraries are working])

	AM_CONDITIONAL([COIN_HAS_ASL],[test a=b])
fi

############
# END PART #
############

AC_LANG_POP(C++)

AC_COIN_VPATH_LINK(src/gmsla_us.run)
AC_COIN_VPATH_LINK(src/gmsdu_us.run)
AC_COIN_VPATH_LINK(src/install/gmsla_un.dat)
AC_COIN_VPATH_LINK(src/install/gmsdu_un.dat)

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/run_tests
                 test/Makefile])

AC_CONFIG_HEADER([src/config_lago.h])

AC_COIN_FINALIZE
