%% Title
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name:RoundVSTruncation.m
%%%% Purpose: To compare Rounding Erro and Truncation Erro
%%%% Aurthor: Yang Yang
%%%% Date: 2015.10.20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control parameters
xspan = [0 10];          % Span of x
K = 18;                 % Times of Change mesh number
N_0 = 8;                % Initial number of span split

%% Compute and statictics erro for first-order derivative

% using first-order back difference format: 
% \frac{f_{i} - f_{i-1}}{\Delta x}
E_11 = zeros(3,K);        % temporary array to store mean absolute erro 
h = zeros(1,K);         
[E_11,h] = FOBDF1st(K,N_0,xspan);  

% using second order centered difference format: 
% \frac{f_{i+1} - f_{i-1}}{2\Delta x}
E_12 = zeros(3,K);        % temporary array to store mean absolute erro        
[E_12,h] = SOCDF1st(K,N_0,xspan);  


%% Compute and statistics erro for second-order derivative

% using three point first-order back difference format: 
% \frac{f_{i} - 2f_{i-1} + f_{i-2}}{\Delta x^{2}}
E_21 = zeros(3,K);        % temporary array to store mean absolute erro          
[E_21,h] = FOBDF2nd(K,N_0,xspan);  

% using three point second order centered difference format: 
% \frac{f_{i+1} - 2f_{i} + f_{i-1}}{\Delta x^{2}}
E_22 = zeros(3,K);        % temporary array to store mean absolute erro          
[E_22,h] = SOCDF2nd(K,N_0,xspan);  


%% plot and compare details
% First order derivative
figureplot1(h, [E_11;E_12])

% Second order derivative
figureplot2(h, [E_21;E_22])

%% Output data
save Result K N_0 xspan h E_11 E_12 E_21 E_22
