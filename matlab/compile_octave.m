% You can also change ../lib for the correspoding install path
% Octave
mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers ...
    --mex --output bayesoptcont.mex bayesoptmex.c

mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers ...
    --mex --output bayesoptdisc.mex bayesoptdiscmex.c

%     if exist('../lib/libbayesopt.a','file')
%     disp('Compiling static library');
%      mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers ...
%         --mex --output bayesoptcont.mex bayesoptmex.c

%     mkoctfile -L../lib -lbayesopt -lnlopt -I../include -I../wrappers ...
%         --mex --output bayesoptdisc.mex bayesoptdiscmex.c
% else % TODO: Does not work in MacOS
%     if (~ismac)
%         disp('Compiling dynamic library');
%         mkoctfile -L../lib -l:bayesopt.so -lnlopt -I../include -I../wrappers ...
%             --mex --output bayesoptcont.mex bayesoptmex.c

%         mkoctfile -L../lib -l:bayesopt.so -lnlopt -I../include -I../wrappers ...
%             --mex --output bayesoptdisc.mex bayesoptdiscmex.c
%     else
%         disp('Dynamic library not supported in MacOS');
%     end
% end
    
