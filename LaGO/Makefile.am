# Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske
# Template from: Andreas Waechter

AUTOMAKE_OPTIONS = foreign

SUBDIRS = src

#if ALWAYS_FALSE
#  SUBDIRS += test
#endif

EXTRA_DIST = \
	README \
	INSTALL \
	LICENSE \
	AUTHORS

.PHONY:

# tmpdir = tmpdir

#for using IpoptCalculatedQuantities
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
  $(IPOPTSRCDIR)/src/Algorithm/IpIpoptCalculatedQuantities.hpp \
  $(IPOPTSRCDIR)/src/Algorithm/IpIpoptNLP.hpp \
  $(IPOPTSRCDIR)/src/Algorithm/IpIpoptData.hpp \
  $(IPOPTSRCDIR)/src/Algorithm/IpNLPScaling.hpp \
  $(IPOPTSRCDIR)/src/Algorithm/IpIteratesVector.hpp \
  $(IPOPTSRCDIR)/src/Algorithm/IpTimingStatistics.hpp \
  $(IPOPTSRCDIR)/src/LinAlg/IpScaledMatrix.hpp \
  $(IPOPTSRCDIR)/src/LinAlg/IpSymScaledMatrix.hpp \
  $(IPOPTSRCDIR)/src/LinAlg/IpCompoundVector.hpp

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

# addlibsfile = lago_addlibs.txt

# install-exec-local:
# 	$(install_sh_DATA) $(addlibsfile) $(DESTDIR)$(libdir)/$(addlibsfile)
# 
# uninstall-local:
# 	rm -f $(DESTDIR)$(libdir)/$(addlibsfile)


# Stuff to clean up

DISTCLEANFILES =
				
include ../BuildTools/Makemain.inc