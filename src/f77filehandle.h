// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef F77FILEHANDLE_H
#define F77FILEHANDLE_H

/** Opens a file for Fortran-routines.
    @param name Filename.
    @return The filedescriptor, you can use in your Fortran-code to access the file. -1, if no free descriptor was available.
*/
int openF77file(char* name);

/** Closes a file, which was opened with openf77file before.
    @param fd The filedescriptor, you got from openf77file.
*/
void closeF77file(int fd);

#endif // F77FILEHANDLE_H
