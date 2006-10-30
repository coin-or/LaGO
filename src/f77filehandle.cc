// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "f77filehandle.h"
#include "standard.h"

#ifdef NOUNDERSCORE
#define F77OPENFILE f77openfile
#define F77CLOSEFILE f77closefile
#else
#define F77OPENFILE f77openfile_
#define F77CLOSEFILE f77closefile_
#endif

/** Opens a file for Fortran.
    @param id The device-number to use.
    @param filename The filename to open.
    @param namelength The length of the filename.
*/
extern "C" void F77OPENFILE(int&, char*, int&);

/** Closes a file for Fortran.
    @param id The device-number of the file to close.
*/
extern "C" void F77CLOSEFILE(int&);

list<int> free_descriptors;
bool initialized=false;

int openF77file(char* name) {
	if (!initialized) {// first run
		for (int i=50; i<90; i++) free_descriptors.push_back(i);
		initialized=true;
	}
	
	if (free_descriptors.empty()) return -1; // no free file descriptor left
		
	int fd=free_descriptors.front();
	free_descriptors.pop_front();
	
	int namelen=strlen(name);
	F77OPENFILE(fd, name, namelen);
	
	return fd;
}

void closeF77file(int fd) {
	F77CLOSEFILE(fd);
	free_descriptors.push_back(fd);
}

