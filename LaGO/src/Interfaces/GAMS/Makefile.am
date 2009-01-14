# Copyright (C) 2007-2009 Stefan Vigerske
# All Rights Reserved.

## $Id$

# Author:  Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

noinst_LTLIBRARIES = libLaGOInterfaceGams.la
libLaGOInterfaceGams_la_SOURCES = LaGOGamsReader.cpp LaGOGamsFunction.cpp LaGOG2DFunction.cpp LaGOGamsSolver.cpp

includecoindir = $(includedir)/coin
includecoin_HEADERS = LaGOGamsFunction.hpp LaGOG2DFunction.hpp LaGOGamsReader.hpp

AM_CPPFLAGS = $(GAMSIO_CPPFLAGS) $(FILIB_CPPFLAGS) $(FORTRAN_NAMEMANGLING_CFLAG) \
 -I`$(CYGPATH_W) $(srcdir)/../../Base` \
 -I`$(CYGPATH_W) $(srcdir)/../../Problem` \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(LAGOSRCDIR)/../Cgc/src` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
 -I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc`
 
if COIN_HAS_COUENNE
AM_CPPFLAGS += -I`$(CYGPATH_W) $(COUENNESRCDIR)/inc` \
  -I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
  -I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
  -I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
  -I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
  -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
  -I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
  -I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
  -I`$(CYGPATH_W) $(CGLOBJDIR)/inc` \
  -I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
  -I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
  -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
  -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
  -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
  -I`$(CYGPATH_W) $(BONMINOBJDIR)/inc`
endif

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../../BuildTools/Makemain.inc
