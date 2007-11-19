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

if ALWAYS_FALSE
  SUBDIRS += test
endif

EXTRA_DIST = \
	README \
	INSTALL \
	LICENSE \
	AUTHORS

# tmpdir = tmpdir

test: all
	cd test; $(MAKE) test

unitTest: test

.PHONY: test unitTest

# Stuff to clean up
clean-local:
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

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