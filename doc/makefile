# Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This code is published under the Common Public License.
#
# Author: Stefan Vigerske

CC  = gcc


doc : paramdochtml doxygen.conf footer.html ../src/*.cc ../src/*.h
	doxygen doxygen.conf
	echo "<HTML><BODY>" > html/paramdoc.html
	bash -c 'for file in ../src/*.h ../src/*.cc ; do paramdochtml < $$file >> html/paramdoc.html; done'
	echo "</BODY></HTML>" >> html/paramdoc.html
	echo "<HTML><BODY>" > html/shortparamdoc.html
	bash -c 'for file in ../src/*.h ../src/*.cc ; do paramdochtml 1 < $$file >> html/shortparamdoc.html; done'
	echo "</BODY></HTML>" >> html/shortparamdoc.html
	echo "<HTML><BODY>" > html/veryshortparamdoc.html
	bash -c 'for file in ../src/*.h ../src/*.cc ; do paramdochtml 2 < $$file >> html/veryshortparamdoc.html; done'
	echo "</BODY></HTML>" >> html/veryshortparamdoc.html

paramdochtml : paramdochtml.lex
	flex -oparamdochtml.c paramdochtml.lex
	$(CC) paramdochtml.c -o paramdochtml

clean:
	rm -f paramdochtml paramdochtml.c
	rm -r html
