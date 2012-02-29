function x_out = langermann(x_in)
% Bounds [3,5]
    a = [3,5,2,1,7];
    b = [5,2,1,4,9];
    c = [1,2,5,2,3];
    
    x = x_in(1); 
    y = x_in(2);
  
    t1 = exp(-(x-a).^2/pi -(y-b).^2/pi);
    t2 = cos(pi*(x-a).^2 + pi*(y-b).^2);
    
    x_out = sum(c.*t1.*t2);