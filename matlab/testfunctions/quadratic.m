function [value] = quadratic(xin)
% Simple quadratic function [0,1]
    xin = xin - 0.53;
    value = dot(xin,xin) + 10.0;
    
return
    