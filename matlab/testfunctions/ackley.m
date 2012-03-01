function y = ackley(x)
% Bounds -32.768, 32.768
% Min x = zeros, y = 0

n = length(x);
a = 20;
b = 0.2;
c = 2*pi;
ccx = cos(c*x);

y = -a*exp(-b*sqrt(sum(x.^2)/n)) - exp(sum(ccx)/n) + a + exp(1);