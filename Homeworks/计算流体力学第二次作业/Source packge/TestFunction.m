function y = TestFunction( x )
% This function is used to compute test function's value
% here input x and y are column
y = zeros(4,length(x));
y(1,:) = x + 1;
y(2,:) = x.^2 + x + 1;
y(3,:) = exp(x);
y(4,:) = sin(5*x) + cos(10*x);

end

