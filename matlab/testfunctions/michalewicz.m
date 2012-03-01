function y = michalewicz(x)
%Bounds [0,pi]
%Min = -4687 (n=5)
    n = length(x);
    m = 1;
    ii = 1:n;
    
    y = -sum(sin(x).*(sin(ii.*x.^2/pi)).^(2*m));