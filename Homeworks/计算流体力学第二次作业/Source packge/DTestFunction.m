function y = DTestFunction( x )
% This function is used to compute test function's first order Derivative value
% here input x and y are column
y = zeros(4,length(x));
y(1,:) = 1;
y(2,:) = 2*x + 1;
y(3,:) = exp(x);
y(4,:) = 5*cos(5*x) - 10*sin(10*x);

end