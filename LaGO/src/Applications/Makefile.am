# Copyright (C) 2007 Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske

EXTRA_DIST = gmsla_us.run install/gmsla_un.dat

AM_CPPFLAGS = \
 $(FILIB_CPPFLAGS) \
 -I`$(CYGPATH_W) $(srcdir)/../Base` \
 -I`$(CYGPATH_W) $(srcdir)/../Problem` \
 -I`$(CYGPATH_W) $(srcdir)/../Algorithms` \
 -I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
 -I`$(CYGPATH_W) $(srcdir)/../Interfaces/GAMS` \
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
 -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
 -I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
 -I`$(CYGPATH_W) $(BONMINOBJDIR)/inc`

if COIN_HAS_CPX
AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` -I`$(CYGPATH_W) $(CPXINCDIR)`
else
AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` -I`$(CYGPATH_W) $(CLPSRCDIR)/src` -I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
endif

lib_LTLIBRARIES = libLaGO.la
libLaGO_la_SOURCES =
libLaGO_la_LIBADD = ../Algorithms/libLaGOAlgorithms.la ../Interfaces/libLaGOInterfaces.la ../Problem/libLaGOProblem.la ../Base/libLaGOBase.la \
 $(RANLIBC_LIB) $(FILIB_LIB)

if COIN_HAS_GAMSIO
bin_PROGRAMS = LaGOGAMS
LaGOGAMS_SOURCES = LaGOGAMS.cpp

LaGOGAMS_LDADD = ../Interfaces/GAMS/libLaGOInterfaceGams.la libLaGO.la

LaGOGAMS_LDADD += \
$(BONMINOBJDIR)/src/CbcBonmin/libbonmin.la \
$(CBCOBJDIR)/src/libCbc.la \
$(CGLOBJDIR)/src/libCgl.la \
$(OSIOBJDIR)/src/libOsi.la \
$(COINUTILSOBJDIR)/src/libCoinUtils.la \
$(IPOPTOBJDIR)/src/Interfaces/libipopt.la \
$(GAMSIO_LIBS) $(SMAG_LIBS) $(GAMSIO_ADDLIBS)

if COIN_HAS_CPX
#to make sure that libOsiCpx comes in the command line before -lcplex, we have to add OsiCpx to ADDLIBS
ADDLIBS += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la
#Cgl seems to insist on OsiClp library, and OsiClp depends on Clp
LaGOGAMS_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la $(CLPOBJDIR)/src/libClp.la
ADDLIBS += `cat $(CLPOBJDIR)/clp_addlibs.txt`
else
LaGOGAMS_LDADD += $(OSIOBJDIR)/src/OsiClp/libOsiClp.la $(CLPOBJDIR)/src/libClp.la
ADDLIBS += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif

ADDLIBS += `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt $(CGLOBJDIR)/cgl_addlibs.txt $(CBCOBJDIR)/cbc_addlibs.txt $(OSIOBJDIR)/osi_addlibs.txt $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt`

LaGOGAMS_LDFLAGS = $(ADDLIBS) -rpath $(libdir)

#lib_LTLIBRARIES = libgamspreprocess.la
#libgamspreprocess_la_SOURCES = gamspreprocess.cc preprocessapi.h
#libgamspreprocess_la_LDFLAGS = -shared -shared-libgcc -Wl,-soname,libgamspreprocess.so -Wl,-export-dynamic
#LaGO_LDFLAGS += -Wl,-export-dynamic

install-exec-local: install-binPROGRAMS
	ln -sf $(DESTDIR)$(bindir)/LaGOGAMS gmsla_ux.out
	zip gmsla_.zip gmsla_ux.out gmsla_us.run install/gmsla_un.dat
	$(install_sh_DATA) gmsla_.zip $(DESTDIR)$(bindir)/gmsla_.zip
	rm -f gmsla_ux.out

uninstall-local:
	rm -f $(DESTDIR)$(bindir)/gmsla_.zip

endif

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../../BuildTools/Makemain.inc
