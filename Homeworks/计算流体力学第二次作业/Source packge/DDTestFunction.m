function y = DDTestFunction( x )
% This function is used to compute test function's second order Derivative value
% here input x and y are column
y = zeros(4,length(x));
y(1,:) = 0;
y(2,:) = 2;
y(3,:) = exp(x);
y(4,:) = -25*sin(5*x) - 100*cos(10*x);

end
