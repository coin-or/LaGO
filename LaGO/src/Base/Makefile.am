# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

CGCSRCDIR = $(LAGOSRCDIR)/../Cgc/src

lib_LTLIBRARIES = libLaGOBase.la
libLaGOBase_la_SOURCES = LaGOSparseVector.cpp LaGODenseVector.cpp LaGOSymSparseMatrix.cpp \
  LaGOFunction.cpp LaGORandomNumber.cpp

includecoindir = $(includedir)/coin
includecoin_HEADERS = ../../inc/config_lago.h LaGOConfig.h LaGObase.hpp LaGOSmartPtr.hpp \
  LaGODenseVector.hpp LaGOSparseVector.hpp \
  LaGOMatrix.hpp LaGOSymMatrix.hpp LaGOSymSparseMatrix.hpp \
  LaGOFunction.hpp LaGOScaledFunction.hpp LaGORandomNumber.hpp \
  $(CGCSRCDIR)/DynNet.h $(CGCSRCDIR)/PLess.h $(CGCSRCDIR)/NodeId.h $(CGCSRCDIR)/NodeIdMgr.h

if COIN_HAS_FILIB
includecoin_HEADERS += LaGOIntervalVector.hpp
libLaGOBase_la_SOURCES += LaGOIntervalVector.cpp
endif

AM_CPPFLAGS = $(RANLIBC_CPPFLAGS) \
 $(FILIB_CPPFLAGS) \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(LAGOSRCDIR)/../Cgc/src` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
 -I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc`
 

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../BuildTools/Makemain.inc
