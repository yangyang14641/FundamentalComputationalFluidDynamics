%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name: FOIM.m
%%%% Program purpose: Analysis Fist-order implicit scheme's accurate
%%%% Aurthor: Yang Yang
%%%% Date: 2015.11.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Fist-order implicit  
% \frac{u_{j}^{n+1}-u_{j}^{n}}{\Delta t}=\frac{u_{j+1}^{n+1}-u_{j-1}^{n+1}}{2\Delta x}

%% basic parameters
% computational parameters
xspan = [0 3];
Nt = 1;
c = 1.1;
% mesh parameters
Nx_0 = 8;
RefineTimes = 11;

MeanErro = zeros(1,RefineTimes+1);
h = zeros(1,RefineTimes+1);

fprintf('Fist-order implicit: c = %f\n',c);
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

% assemble matrix
A = zeros(Nx+1,Nx+1);

A(1,Nx) = -1/2*c;
A(1,1) = 1;
A(1,2) = 1/2*c;
for i = 2:Nx
    A(i,i-1) = -1/2*c;
    A(i,i) = 1;
    A(i,i+1) = 1/2*c;
end
A(Nx+1,Nx) = -1/2*c;
A(Nx+1,Nx+1) = 1;
A(Nx+1,2) = 1/2*c;


% forward in time step
for n = 2:Nt+1
    u(n,:) = (A\(u(n-1,:)'))';    % All points using Matrix Algebra
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

    % assemble matrix
    A = zeros(Nx+1,Nx+1);

    A(1,Nx) = -1/2*c;
    A(1,1) = 1;
    A(1,2) = 1/2*c;
    for i = 2:Nx
        A(i,i-1) = -1/2*c;
        A(i,i) = 1;
        A(i,i+1) = 1/2*c;
    end
    A(Nx+1,Nx) = -1/2*c;
    A(Nx+1,Nx+1) = 1;
    A(Nx+1,2) = 1/2*c;


    % forward in time step
    for n = 2:Nt+1
        u(n,:) = (A\(u(n-1,:)'))';    % All points using Matrix Algebra
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
title('First-order implicit');