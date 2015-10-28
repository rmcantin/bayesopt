# -------------------------------------------------------------------------
#    This file is part of BayesOpt, an efficient C++ library for 
#    Bayesian optimization.
#
#    Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
#
#    BayesOpt is free software: you can redistribute it and/or modify it 
#    under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    BayesOpt is distributed in the hope that it will be useful, but 
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------


MATLABROOT= C:\Program Files (x86)\MATLAB\R2010a
CC=mingw32-c++
INCLUDES=-I"$(MATLABROOT)\extern\include" -I"." -I"..\include" -I"..\wrappers"
LIBS=-L"$(MATLABROOT)\bin\win32" -lmex -lmx -leng -lmat
FLAGS=-shared -DMATLAB_MEX_FILE -Wl,--export-all-symbols
BOPT=libbayesopt.dll

all: disc cont

cont:
	$(CC) $(FLAGS) $(INCLUDES) -o bayesoptcont.mexw32  bayesoptmex.c $(BOPT) $(LIBS)

disc:
	$(CC) $(FLAGS) $(INCLUDES) -o bayesoptdisc.mexw32  bayesoptdiscmex.c $(BOPT) $(LIBS)










