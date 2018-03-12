function [x,delta] = meshfun2( x_span, N )
% This function is designed for generate property mesh
factor = 1.1;                              % space factor

delta = zeros(1,N-1);                      % Temporaray array to store delta
x = zeros(1,N);                            % store mesh node
a1 = (x_span(2) - x_span(1)) ...
    * (1 - factor) / (1- factor^(N-1));    % compute fist element in delta

delta(1) = a1;                             %Initial first element in delta array
x(1) = x_span(1);                          %Initial fist element in x array

for i = 2:N-1                              % generate equlibrium distribution nodes
     delta(i) = delta(i-1)*factor;
     x(i) = x(i-1) + delta(i-1);
end
x(i+1) = x(i) + delta(i);                  % The last element in mesh nodes

end