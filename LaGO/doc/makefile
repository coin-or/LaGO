# Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This code is published under the Common Public License.
#
# Author: Stefan Vigerske

CC  = gcc

doc : paramdochtml ../src/*.cc ../src/*.h
	echo "<HTML><BODY>" > paramdoc.html
	bash -c 'for file in ../src/*.h ../src/*.cc ; do paramdochtml < $$file >> paramdoc.html; done'
	echo "</BODY></HTML>" >> paramdoc.html
	echo "<HTML><BODY>" > shortparamdoc.html
	bash -c 'for file in ../src/*.h ../src/*.cc ; do paramdochtml 1 < $$file >> shortparamdoc.html; done'
	echo "</BODY></HTML>" >> shortparamdoc.html
	echo "<HTML><BODY>" > veryshortparamdoc.html
	bash -c 'for file in ../src/*.h ../src/*.cc ; do paramdochtml 2 < $$file >> veryshortparamdoc.html; done'
	echo "</BODY></HTML>" >> veryshortparamdoc.html

paramdochtml : paramdochtml.lex
	flex -oparamdochtml.c paramdochtml.lex
	$(CC) paramdochtml.c -o paramdochtml

clean:
	rm -f paramdochtml paramdochtml.c paramdoc.html shortparamdoc.html veryshortparamdoc.html
