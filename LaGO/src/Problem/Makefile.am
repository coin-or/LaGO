# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libLaGOProblem.la
libLaGOProblem_la_SOURCES = LaGOMINLPData.cpp LaGOBlockFunction.cpp LaGOSparsity.cpp LaGOSampleSet.cpp LaGOBoxMinimizationProblem.cpp

include_HEADERS = LaGOMINLPData.hpp LaGOBlockFunction.hpp LaGOSparsity.hpp LaGOCurvature.hpp LaGOQuadraticFunction.hpp LaGOSampleSet.hpp

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
 -I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc`
 

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../BuildTools/Makemain.inc
