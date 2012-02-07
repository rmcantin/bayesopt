function [value] = testfunc(xin)
% Test function
    disp(size(xin));
    xin = xin - 0.53;
    value = dot(xin,xin) + 10.0;
    
return
    