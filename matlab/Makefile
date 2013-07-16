MATLABROOT=C:\Program Files (x86)\MATLAB\R2010a
CC=mingw32-c++
INCLUDES=-I"$(MATLABROOT)\extern\include" -I"." -I"..\include" -I"..\wrappers"
LIBS=-L"$(MATLABROOT)\bin\win32" -lmex -lmx -leng -lmat
FLAGS=-shared -DMATLAB_MEX_FILE -Wl,--export-all-symbols

all: disc cont

cont:
	$(CC) $(FLAGS) $(INCLUDES) -o bayesoptcont.mexw32  bayesoptmex.c libbayesopt.dll $(LIBS)

disc:
	$(CC) $(FLAGS) $(INCLUDES) -o bayesoptdisc.mexw32  bayesoptdiscmex.c libbayesopt.dll $(LIBS)

#mingw32-c++ -shared -I"C:\Program Files (x86)\MATLAB\R2010a\extern\include" -I"." -I"..\include" -I"..\wrappers" -DMATLAB_MEX_FILE -o bayesoptdisc.mexw32 -Wl,--export-all-symbols bayesoptdiscmex.c libbayesopt.dll -L"C:\Program Files (x86)\MATLAB\R2010a\bin\win32" -lmex -lmx -leng -lmat









