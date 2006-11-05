# Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske
# Template from: Andreas Waechter

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST =

bin_PROGRAMS = LaGO
LaGO_SOURCES = bcp.cc bcp.h boxfind.cc boxfind.h column.cc column.h config_lago.h cuts.cc cuts.h decomp.cc decomp.h func.cc func.h graph.h ipopt2.cc ipopt2.h lagheu.cc lagheu.h LaGO.cc LaGOConfig.h linrelax.cc linrelax.h MINLP.cc MINLP.h MINLPData.cc MINLPData.h minlpopt.cc minlpopt.h MINLPView.cc MINLPView.h node.cc node.h opt.cc opt.h osi.cc osi.h param.cc param.h polynom.cc polynom.h problem.cc problem.h quaduest.cc quaduest.h relax.cc relax.h relaxopt.cc relaxopt.h rmp.cc rmp.h sampling.cc sampling.h standard.h tools.cc tools.h usermatrix.cc usermatrix.h uservector.h \
gams.cc gams.h preprocessapi.h ampl.cc ampl.h

LaGO_LDADD = $(COINUTILSOBJDIR)/src/libCoinUtils.la \
$(CLPOBJDIR)/src/libClp.la \
$(CGLOBJDIR)/src/libCgl.la \
$(OSIOBJDIR)/src/libOsi.la \
$(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
$(IPOPTOBJDIR)/src/Interfaces/libipopt.la
#add OsiCpx if Cplex available
#$(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la

LIBS += $(ADDLIBS) `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt $(CLPOBJDIR)/clp_addlibs.txt $(CGLOBJDIR)/cgl_addlibs.txt $(OSIOBJDIR)/osi_addlibs.txt $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt`

LaGO_LDFLAGS = -rpath $(libdir)

AM_CPPFLAGS = \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
 -I`$(CYGPATH_W) $(COINUTILSSRCDIR)/inc` \
 -I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
 -I`$(CYGPATH_W) $(CGLSRCDIR)/inc` \
 -I`$(CYGPATH_W) $(OSISRCDIR)/src` \
 -I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
 -I`$(CYGPATH_W) $(OSISRCDIR)/inc` \
 -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
 -I`$(CYGPATH_W) $(CLPSRCDIR)/inc` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
 -I`$(CYGPATH_W) $(IPOPTSRCDIR)/inc`

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES =

include ../../BuildTools/Makemain.inc
