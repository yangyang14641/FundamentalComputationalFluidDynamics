%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name: LeapFrog.m
%%%% Program purpose: Analysis Leap-frog scheme's accurate
%%%% Aurthor: Yang Yang
%%%% Date: 2015.11.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \frac{u_{j}^{n+1}-u_{j}^{n-1}}{2\Delta t}=\frac{u_{j+1}^{n}-u_{j-1}^{n}}{2\Delta x}

%% basic parameters
% computational parameters
xspan = [0 3];
Nt = 1;
c = 0.5;
% mesh parameters
Nx_0 = 8;
RefineTimes = 11;

MeanErro = zeros(1,RefineTimes+1);
h = zeros(1,RefineTimes+1);

fprintf('LeapFrog: c = %f\n',c);
%% Compute solution in Initial mesh 
Nx = Nx_0;
h(1) = (xspan(2) - xspan(1)) / Nx;
fprintf('Current mesh nodes amount: %d\n',Nx)

% initial variables
delta_x = (xspan(2) - xspan(1))/Nx;
delta_t = c*delta_x;
tspan = [0, Nt*delta_t];

x = zeros(1,Nx+1);
t = zeros(1,Nt+1);
u = zeros(Nt+1,Nx+1);
u_real = zeros(1,Nx+1);

% initial mesh nodes
x = xspan(1):delta_x:xspan(2);   % x(i) spatial nodes
t = tspan(1):delta_t:tspan(2);   % t(i) time nodes

u(1,:) = sin(2*pi*(x));          % u(i,0) initial u using initial condition
u_temp = sin(2*pi*(x+delta_t));  % u(i,0) visual initial condition

% forward in time step
n = 2;
for i = 2:Nx
    u(n,i) = u_temp(i) - c*( u(n-1,i+1) - u(n-1,i-1) );        % inner points
end
u(n,1) = u_temp(1) - c*( u(n-1,2) - u(n-1,Nx) );             % boundary points
u(n,Nx+1) = u_temp(Nx+1) - c*( u(n-1,2) - u(n-1,Nx) );    % boundary points


for n = 3:Nt+1
    for i = 2:Nx
        u(n,i) = u(n-2,i) - c*( u(n-1,i+1) - u(n-1,i-1) );        % inner points
    end
    u(n,1) = u(n-2,1) - c*( u(n-1,2) - u(n-1,Nx) );             % boundary points
    u(n,Nx+1) = u(n-2,Nx+1) - c*( u(n-1,2) - u(n-1,Nx) );    % boundary points
end

u_real(:) = sin(2*pi*(x - delta_t*Nt));
MeanErro(1) = norm(u(Nt+1,:)-u_real,1);

%% Compute solution and refine mesh
parfor k = 1:RefineTimes
    xspan = [0 3];
    Nt = 1;
    N_step = 2^(k);
    Nx = Nx_0*N_step;
    h(k+1) = (xspan(2) - xspan(1)) / Nx;
    fprintf('Current mesh nodes amount: %d\n',Nx)

    % initial variables
    delta_x = (xspan(2) - xspan(1))/Nx;
    delta_t = c*delta_x;
    tspan = [0, Nt*delta_t];


    x = zeros(1,Nx+1);
    t = zeros(1,Nt+1);
    u = zeros(Nt+1,Nx+1);
    u_real = zeros(1,Nx+1);
    
    % initial mesh nodes
    x = xspan(1):delta_x:xspan(2);   % x(i) spatial nodes
    t = tspan(1):delta_t:tspan(2);   % t(i) time nodes

    u(1,:) = sin(2*pi*(x));          % u(i,0) initial u using initial condition
    u_temp = sin(2*pi*(x+delta_t));  % u(i,0) visual initial condition

% forward in time step
n = 2;
for i = 2:Nx
    u(n,i) = u_temp(i) - c*( u(n-1,i+1) - u(n-1,i-1) );        % inner points
end
u(n,1) = u_temp(1) - c*( u(n-1,2) - u(n-1,Nx) );             % boundary points
u(n,Nx+1) = u_temp(Nx+1) - c*( u(n-1,2) - u(n-1,Nx) );    % boundary points


for n = 3:Nt+1
    for i = 2:Nx
        u(n,i) = u(n-2,i) - c*( u(n-1,i+1) - u(n-1,i-1) );        % inner points
    end
    u(n,1) = u(n-2,1) - c*( u(n-1,2) - u(n-1,Nx) );             % boundary points
    u(n,Nx+1) = u(n-2,Nx+1) - c*( u(n-1,2) - u(n-1,Nx) );    % boundary points
end

      u_real(:) = sin(2*pi*(x - delta_t*Nt));
      MeanErro(k+1) = norm(u(Nt+1,:)-u_real,1);
end

%% release memory
clear u u_real;

%% plot
figure('Color',[1 1 1]);
loglog(h,MeanErro,'linewidth',2);
grid on;
xlabel('lg(h)');
ylabel('lg(E)');
title('LeapFrog');