function y = NACA0012( x , x0)
% This function generate NACA0012 Airfoil Shape

x = x - x0;
y = 0.6*(0.2969*sqrt(x) - 0.126*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4);

end

