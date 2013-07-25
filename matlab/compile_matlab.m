% You can also change ../lib for the correspoding install path
% MATLAB
if (ispc)
    if exist('../bin/Release/bayesopt.dll','file')
        disp('Compiling dynamic library');
        mex -DBAYESOPT_DLL -output bayesoptcont bayesoptmex.c ...
            -L..\lib\Release -L. -lbayesopt ...
            -I..\include -I..\wrappers
        mex -DBAYESOPT_DLL -output bayesoptdisc bayesoptdiscmex.c ...
            -L..\lib\Release -L. -lbayesopt ...
            -I..\include -I..\wrappers
    else
        disp('Compiling static library');
        mex -output bayesoptcont bayesoptmex.c ...
            -L../lib/Release -lbayesopt -lnlopt ...
            -I../include -I../wrappers
        
        mex -output bayesoptdisc bayesoptdiscmex.c ...
            -L../lib/Release -lbayesopt -lnlopt ...
            -I../include -I../wrappers
    end
else
    mex -output bayesoptcont bayesoptmex.c -L../lib -lbayesopt ...
        -lnlopt -I../include -I../wrappers -I../nlopt/api 

    mex -output bayesoptdisc bayesoptdiscmex.c -L../lib -lbayesopt ...
        -lnlopt -I../include -I../wrappers -I../nlopt/api 

    % if exist('../lib/libbayesopt.a','file')
    %     disp('Compiling static library');
    %     mex -output bayesoptcont bayesoptmex.c ../lib/libbayesopt.a ...
    %     ../lib/libnlopt.a -I../include -I../wrappers -I../nlopt/api 

    %     mex -output bayesoptdisc bayesoptdiscmex.c ../lib/libbayesopt.a ...
    %         ../lib/libnlopt.a -I../include -I../wrappers -I../nlopt/api 
    % else
    %     if exist('../lib/bayesopt.so','file')
    %         disp('Compiling dynamic library');
    %         mex -g -output bayesoptcont bayesoptmex.c ../lib/bayesopt.so ...
    %             -I../include -I../wrappers

    %         mex -g -output bayesoptdisc bayesoptdiscmex.c ../lib/bayesopt.so ...
    %             -I../include -I../wrappers
                
    %     else
    %         disp('Error: File not found');
    %     end
    % end
end