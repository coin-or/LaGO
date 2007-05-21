# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske

#EXTRA_DIST = gmsla_us.run install/gmsla_un.dat
EXTRA_DIST = 

bin_PROGRAMS = LaGO
LaGO_SOURCES = 

AM_CPPFLAGS = \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
 -I`$(CYGPATH_W) $(CGLOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
 -I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
 -I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc`

LaGO_LDADD =
 
if COIN_HAS_CPX
AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` -I`$(CYGPATH_W) $(CPXINCDIR)`
LaGO_LDADD += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
#Cgl seem to insist on OsiClp library, and OsiClp depends on Clp
LaGO_LDADD += $(CLPOBJDIR)/src/libClp.la $(OSIOBJDIR)/src/OsiClp/libOsiClp.la
LIBS += `cat $(CLPOBJDIR)/clp_addlibs.txt`
else
AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` -I`$(CYGPATH_W) $(CLPSRCDIR)/src` -I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
LaGO_LDADD += $(CLPOBJDIR)/src/libClp.la $(OSIOBJDIR)/src/OsiClp/libOsiClp.la
LIBS += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif

LaGO_LDADD += \
$(OSIOBJDIR)/src/libOsi.la \
$(CGLOBJDIR)/src/libCgl.la \
$(COINUTILSOBJDIR)/src/libCoinUtils.la \
$(IPOPTOBJDIR)/src/Interfaces/libipopt.la

ADDLIBS += `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt $(CGLOBJDIR)/cgl_addlibs.txt $(OSIOBJDIR)/osi_addlibs.txt $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt`

LaGO_LDFLAGS = $(ADDLIBS) -rpath $(libdir) 

if COIN_HAS_GAMSIO
#LaGO_SOURCES += gams.cc gams.h preprocessapi.h 

#lib_LTLIBRARIES = libgamspreprocess.la
#libgamspreprocess_la_SOURCES = gamspreprocess.cc preprocessapi.h
#libgamspreprocess_la_LDFLAGS = -shared -shared-libgcc -Wl,-soname,libgamspreprocess.so -Wl,-export-dynamic
#LaGO_LDFLAGS += -Wl,-export-dynamic

#install-exec-local: install-binPROGRAMS
#	ln -sf $(DESTDIR)$(bindir)/LaGO gmsla_ux.out
#	zip gmsla_.zip gmsla_ux.out gmsla_us.run gmsdu_us.run install/gmsla_un.dat install/gmsdu_un.dat
#	$(install_sh_DATA) gmsla_.zip $(DESTDIR)$(bindir)/gmsla_.zip
#	rm -f gmsla_ux.out

#uninstall-local:
#	rm -f $(DESTDIR)$(bindir)/gmsla_.zip

else
#LaGO_SOURCES += ampl.cc ampl.h
endif

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../BuildTools/Makemain.inc
