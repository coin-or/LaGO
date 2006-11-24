# Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske
# Template from: Andreas Waechter

AUTOMAKE_OPTIONS = foreign

SUBDIRS = interval

# if ALWAYS_FALSE
#   SUBDIRS += test
# endif

EXTRA_DIST = \
	get.FILIB++ \
	filib.patch filib.patch2 filib.patch3

.PHONY:

# tmpdir = tmpdir

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

# install-exec-local:
# 	$(install_sh_DATA) $(addlibsfile) $(DESTDIR)$(libdir)/$(addlibsfile)
# 
# uninstall-local:
# 	rm -f $(DESTDIR)$(libdir)/$(addlibsfile)


# Stuff to clean up

DISTCLEANFILES =
				
include ../../BuildTools/Makemain.inc
