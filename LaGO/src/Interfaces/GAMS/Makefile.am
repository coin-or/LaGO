# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.

## $Id: Makefile.am 18 2007-04-02 16:25:41Z stefan $

# Author:  Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libLaGOInterfaceGams.la
libLaGOInterfaceGams_la_SOURCES = LaGOGamsReader.cpp LaGOGamsFunction.cpp

include_HEADERS = LaGOGamsFunction.hpp LaGOGamsReader.hpp

AM_CPPFLAGS = $(GAMSIO_CPPFLAGS) \
 -I`$(CYGPATH_W) $(srcdir)/../../Base` \
 -I`$(CYGPATH_W) $(srcdir)/../../Problem` \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
 -I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc`
 
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../../BuildTools/Makemain.inc
