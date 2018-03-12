function r = meshfun(r_span,N_r)
% This function control mesh nodes distribution
factor = 1.02;
r = zeros(1,N_r);         % array to store r;
r(1) = r_span(1);         % Initial value for r
delta = zeros(1,N_r-1);   % array to store delta;
delta(1) = (r_span(2) - r_span(1))*(1 - factor)/(1 - factor^(N_r-1));

for i = 2:N_r -1
     delta(i) = factor*delta(i-1);
     r(i) = r(i-1) + delta(i-1);
end
r(N_r) = r(N_r-1) + delta(N_r-1); 

end

