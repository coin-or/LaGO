# Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
# All Rights Reserved.
# This code is published under the Common Public License.
#
# Author: Stefan Vigerske

CPP = g++

CPPFLAGS    = -DNO_NULLPOINTER_CHECK -DTNT_NO_BOUNDS_CHECK -DNO_SPARSEMATRIX_ASSERTS -DNDEBUG -O3
#CPPFLAGS	= -g

################ LaGO objects
OBJ		= tools.o func.o usermatrix.o problem.o opt.o param.o MINLPData.o MINLPView.o MINLP.o relaxopt.o relax.o linrelax.o cuts.o boxfind.o decomp.o polynom.o sampling.o minlpopt.o node.o bcp.o column.o rmp.o lagheu.o
#subgradopt.o dualqqp.o 

################ Third Party codes
INC		+= -I../ThirdParty/TNT

INC		+= -I../ThirdParty/ranlib.c/src
LIBS	+= -L../ThirdParty/ranlib.c -lranlib

CPPFLAGS	+= -DMETIS_AVAILABLE
INC		+= -I../ThirdParty/metis-4.0/Lib
LIBS	+= -L../ThirdParty/metis-4.0 -lmetis

CPPFLAGS	+= -DAMPL_AVAILABLE
INC		+= -I../ThirdParty/solvers
LIBS	+= -ldl ../ThirdParty/solvers/amplsolver.a
OBJ		+= ampl.o

CPPFLAGS	+= -DCOIN_AVAILABLE
COINDIR = /home/stefan/soft/COIN/debug
INC		+= -I$(COINDIR)/include
LIBS	+= -L$(COINDIR)/lib -lCgl -lOsiClp -lClp -lOsi -lCoinUtils
OBJ		+= osi.o

CPPFLAGS	+= -DIPOPT_AVAILABLE
LIBS	+= -lgfortranbegin -lgfortran -lg2c -lipopt
OBJ		+= ipopt2.o

#CPPFLAGS	+= -DFILIB_AVAILABLE
#FILIBDIR =
#LIBS	+= -L$(FILIBDIR)/lib -lprim
#INC		+= -I$(FILIBDIR)/include

#CPPFLAGS	+= -DSNOPT_AVAILABLE
#SNOPTDIR =
#LIBS	+= -L$(SNOPTDIR) -lsnopt
#OBJ		+= snopt.o f77files.o f77filehandle.o -lg2c

################ targets

.cc.o:
	$(CPP) $(CPPFLAGS) $(INC) -c $*.cc

LaGO : $(OBJ) LaGO.o
	$(CPP) -o LaGO $^ $(LIBS)

dep depfile:
	$(CPP) -MM $(CPPFLAGS) $(INC) *.cc >depfile

include depfile

clean:
	rm -f depfile *.o LaGO

