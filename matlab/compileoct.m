% You can also change ../lib for the correspoding install path
% Octave
mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers -I../nlopt/api --mex "-Wl,-rpath=../lib" --output bayesopt.mex  bayesoptmex.c
    
