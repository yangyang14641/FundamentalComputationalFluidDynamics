function [E,h] = SOCDF1st(K,N_0,xspan)
% This function is to compute Erro of each mesh density
%% Initial Output
E = zeros(3,K); 
h = zeros(1,K); 
%% Span split
N = N_0*2^(K-1);                           % number of span split 
delta_x = (xspan(2) - xspan(1)) / N;       % step size of x

%% nodes' position
x = zeros(1,N+1);                          % Store noeds' position
x(1) = xspan(1);
for i = 2:N+1
   x(i) = x(i-1) + delta_x;                % Nodes' position
end
clear i                                    % free memory

%% Compute function value on this points
f = TestFunction( x ); 
df = DTestFunction( x );
%ddf = DDTestFunction( x );

%% Main loop to compute mean absolute erro
for k = 1:4
     step = 0;
  for i = 1:K                    % index i represent mesh density level
     temp = 0;
     % normal case
     index = 2^(step);           % define index, using index to replace 1 at each step
     delta = delta_x*index;      % mesh step size in current step
     counter = 0;                % to record how many points we compute
     for j = 1:index:N           % index j represent x mesh nodes
        if j >1 && j<N
             temp = temp + abs( (f(k,j+index) - f(k,j-index)) / (2*delta) - df(k,j) );   % SOCDF
             counter = counter + 1;
        end
     end
     % store
     E(k,i) = temp/counter;
     h(i) = delta;
     step = step + 1;
  end
  
end


end



