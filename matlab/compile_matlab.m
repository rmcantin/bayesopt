% You can also change ../lib for the correspoding install path
% MATLAB
if (ispc)
    mex -output bayesopt bayesoptmex.c ..\lib\Release\bayesopt.lib ...
        ..\lib\Release\nlopt.lib ...
        -I..\include -I..\wrappers -I..\nlopt\api 

   mex -output bayesoptdisc bayesoptdiscmex.c ../lib/Release/bayesopt.lib ...
        ../lib/Release/nlopt.lib -I../include -I../wrappers -I../nlopt/api 
else
    if exist('../lib/libbayesopt.a','file')
        disp('Compiling static library');
        mex -output bayesopt bayesoptmex.c ../lib/libbayesopt.a ...
        ../lib/libnlopt.a -I../include -I../wrappers -I../nlopt/api 

        mex -output bayesoptdisc bayesoptdiscmex.c ../lib/libbayesopt.a ...
            ../lib/libnlopt.a -I../include -I../wrappers -I../nlopt/api 
    else
        if exist('../lib/bayesopt.so','file')
            disp('Compiling dynamic library');
            mex -g -output bayesopt bayesoptmex.c ../lib/bayesopt.so ...
                -I../include -I../wrappers

            mex -g -output bayesoptdisc bayesoptdiscmex.c ../lib/bayesopt.so ...
                -I../include -I../wrappers
                
        else
            disp('Error: File not found');
        end
    end
end