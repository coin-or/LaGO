# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

lib_LTLIBRARIES = libranlib.c.la
libranlib_c_la_SOURCES = ranlib.c/src/ranlib.c ranlib.c/src/com.c ranlib.c/linpack/linpack.c
include_HEADERS = ranlib.c/src/ranlib.h

libranlib_c_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/ranlib.c/src` 

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../BuildTools/Makemain.inc
