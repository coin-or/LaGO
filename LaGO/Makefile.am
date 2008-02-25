# Copyright (C) 2006-2008 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id$
#
# Author: Stefan Vigerske

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

if COIN_HAS_GAMSIO
if COIN_HAS_GAMSSYSTEM
gams-install: install
	echo "Installation of LaGO in GAMS distribution at @GAMS_PATH@."
	cp $(DESTDIR)$(bindir)/gmsla_.zip @GAMS_PATH@
	cd @GAMS_PATH@; ./gamsinst -a
else
gams-install: 
	echo "Installation of LaGO in GAMS distribution skipped. No GAMS system available."
	exit 0
endif
else
gams-install: 
	echo "Installation of LaGO in GAMS distribution skipped. LaGO was build with AMPL interface."
	exit 1
endif

.PHONY: test unitTest gams-install

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