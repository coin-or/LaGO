# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libLaGOProblem.la
libLaGOProblem_la_SOURCES = LaGOMINLPData.cpp LaGOSparsity.cpp LaGOSampleSet.cpp \
  LaGOBlockFunction.cpp LaGOSimpleBlockFunction.cpp LaGOSumFunction.cpp LaGORestrictedFunction.cpp LaGOQuadraticFunction.cpp \
  LaGOBoxMinimizationProblem.cpp LaGOQuadraticOrConvexApproximation.cpp

include_HEADERS = LaGOMINLPData.hpp LaGOSparsity.hpp LaGOCurvature.hpp \
  LaGOBlockFunction.hpp LaGOQuadraticFunction.hpp LaGOQuadraticEstimator.hpp \
  LaGOSampleSet.hpp

AM_CPPFLAGS = \
 $(FILIB_CPPFLAGS) \
 -I`$(CYGPATH_W) $(srcdir)/../Base` \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
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

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../BuildTools/Makemain.inc
