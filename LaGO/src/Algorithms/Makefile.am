# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libLaGOAlgorithms.la
libLaGOAlgorithms_la_SOURCES = \
	LaGOAlgorithm.cpp \
	LaGODecomposition.cpp \
	LaGOSampling.cpp LaGOCurvatureCheck.cpp \
	LaGOBoxReduction.cpp LaGOConstraintPropagation.cpp \
	LaGOQuadraticEstimation.cpp LaGOConvexification.cpp

include_HEADERS = 

AM_CPPFLAGS = $(FILIB_CPPFLAGS) \
 -I`$(CYGPATH_W) $(srcdir)/../Base` \
 -I`$(CYGPATH_W) $(srcdir)/../Problem` \
 -I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
 -I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(LAGOSRCDIR)/../Cgc/src` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
 -I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
 -I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
 -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
 -I`$(CYGPATH_W) $(BONMINOBJDIR)/inc`

if COIN_HAS_CPX
AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
 -I`$(CYGPATH_W) $(CPXINCDIR)`
else
AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
 -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
 -I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
endif

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../BuildTools/Makemain.inc
