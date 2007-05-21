# Copyright (C) Stefan Vigerske 2007
# All Rights Reserved.

## $Id: Makefile.am 18 2007-04-02 16:25:41Z stefan $

# Author:  Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libLaGOBase.la
libLaGOBase_la_SOURCES =  

include_HEADERS = ../../inc/config_lago.h LaGOConfig.h base.hpp

AM_CPPFLAGS = \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`
 

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../BuildTools/Makemain.inc
