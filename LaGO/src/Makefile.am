# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

SUBDIRS = Base Problem

if COIN_HAS_GAMSIO
SUBDIRS += Interfaces/GAMS
endif

SUBDIRS += Applications

EXTRA_DIST = 

#DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../BuildTools/Makemain.inc
