%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name: OneDimensionWave Equation
%%%% Program purpose: Compute and Analysis different numerical scheme 
%%%% Aurthor: Yang Yang
%%%% Date: 2015.11.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fist-order upwind  
%\frac{u_{j}^{n+1}-u_{j}^{n}}{\Delta t}=\frac{u_{j}^{n}-u_{j-1}^{n}}{\Delta x}
% basic parameters
xspan = [0 3];
tspan = [0 3];
Nx = 100;
Nt = 80;
plots = 5;

% initial variables
delta_x = (xspan(2) - xspan(1))/Nx;
delta_t = (tspan(2) - tspan(1))/Nt;
c = delta_t / delta_x;

fprintf('First-order upwind: c = %f\n',c);
x = zeros(1,Nx+1);
t = zeros(1,Nt+1);
u = zeros(Nt+1,Nx+1);

% initial mesh nodes
x = xspan(1):delta_x:xspan(2);   % x(i) spatial nodes
t = tspan(1):delta_t:tspan(2);   % t(i) time nodes

u(1,:) = sin(2*pi*(x));          % u(i,0) initial u using initial condition

% forward in time step
for n = 2:Nt+1
    for i = 2:Nx
        u(n,i) = u(n-1,i) - c*( u(n-1,i) - u(n-1,i-1) );        % inner points
    end
    u(n,1) = u(n-1,1) - c*( u(n-1,1) - u(n-1,Nx) );             % boundary points
    u(n,Nx+1) = u(n-1,Nx+1) - c*( u(n-1,Nx+1) - u(n-1,Nx) );    % boundary points
end
    
% plot results
% surf
f1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',f1);
colormap('cool');
[X,T] = meshgrid(x,t);
surf(T,X,u,'Parent',axes1,'EdgeAlpha',0.8,'EdgeColor','none');
box on;
grid on;
xlabel('t');
ylabel('x');
zlabel('u');
set(axes1,'FontSize',14,'LineWidth',2);
title('Fist-order upwind');

% curves
figure('Color',[1 1 1]);
colormap jet

Nt_temp = 1;
name_cell = cell(1,plots);
for i = 1:plots
    delta_Nt = floor((Nt+1) / plots);
    plot(x,u(Nt_temp,:));
    hold on
    name_cell{i} = ['t = ',num2str(t(Nt_temp))];
    Nt_temp = Nt_temp + delta_Nt;
end
hold off
legend(name_cell);
xlabel('x');
ylabel('u');    
title('Fist-order upwind')

%% Implicit Scheme 
% \frac{u_{j}^{n+1}-u_{j}^{n}}{\Delta t}=\frac{u_{j+1}^{n+1}-u_{j-1}^{n+1}}{2\Delta x}

% basic parameters
xspan = [0 3];
tspan = [0 3];
Nx = 100;
Nt = 80;
plots = 5;

% initial variables
delta_x = (xspan(2) - xspan(1))/Nx;
delta_t = (tspan(2) - tspan(1))/Nt;
c = delta_t / delta_x;

fprintf('Implicit Scheme: c = %f\n',c);
x = zeros(1,Nx+1);
t = zeros(1,Nt+1);
u = zeros(Nt+1,Nx+1);

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

% plot results
% surf
f1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',f1);
colormap('cool');
[X,T] = meshgrid(x,t);
surf(T,X,u,'Parent',axes1,'EdgeAlpha',0.8,'EdgeColor','none');
box on;
grid on;
xlabel('t');
ylabel('x');
zlabel('u');
set(axes1,'FontSize',14,'LineWidth',2);
title('Fist-order implicit')

% curves
figure('Color',[1 1 1]);
colormap jet

Nt_temp = 1;
name_cell = cell(1,plots);
for i = 1:plots
    delta_Nt = floor((Nt+1) / plots);
    plot(x,u(Nt_temp,:));
    hold on
    name_cell{i} = ['t = ',num2str(t(Nt_temp))];
    Nt_temp = Nt_temp + delta_Nt;
end
hold off
legend(name_cell);
xlabel('x');
ylabel('u'); 
title('Fist-order implicit')

%% Lax Scheme
% \frac{u_{j}^{n+1}-\frac{1}{2}\left(  u_{j+1}^{n}+u_{j-1}^{n}\right)  }{\Delta t}=\frac{u_{j+1}^{n}-u_{j-1}^{n}}{2\Delta x}%
% basic parameters
xspan = [0 3];
tspan = [0 3];
Nx = 100;
Nt = 90;
plots = 5;

% initial variables
delta_x = (xspan(2) - xspan(1))/Nx;
delta_t = (tspan(2) - tspan(1))/Nt;
c = delta_t / delta_x;

fprintf('Lax Scheme: c = %f\n',c);
x = zeros(1,Nx+1);
t = zeros(1,Nt+1);
u = zeros(Nt+1,Nx+1);

% initial mesh nodes
x = xspan(1):delta_x:xspan(2);   % x(i) spatial nodes
t = tspan(1):delta_t:tspan(2);   % t(i) time nodes

u(1,:) = sin(2*pi*(x));          % u(i,0) initial u using initial condition

% forward in time step
for n = 2:Nt+1
    for i = 2:Nx
        u(n,i) = 0.5*(1-c)*u(n-1,i+1) - 0.5*(1+c)*u(n-1,i-1);      % inner points
    end
    u(n,1) = 0.5*(1-c)*u(n-1,2) - 0.5*(1+c)*u(n-1,Nx);                % boundary points
    u(n,Nx+1) = 0.5*(1-c)*u(n-1,2) - 0.5*(1+c)*u(n-1,Nx);       % boundary points
end
    
% plot results
% surf
f1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',f1);
colormap('cool');
[X,T] = meshgrid(x,t);
surf(T,X,u,'Parent',axes1,'EdgeAlpha',0.8,'EdgeColor','none');
box on;
grid on;
xlabel('t');
ylabel('x');
zlabel('u');
set(axes1,'FontSize',14,'LineWidth',2);
title('Lax Scheme')

% curves
figure('Color',[1 1 1]);
colormap jet

Nt_temp = 1;
name_cell = cell(1,plots);
for i = 1:plots
    delta_Nt = floor((Nt+1) / plots);
    plot(x,u(Nt_temp,:));
    hold on
    name_cell{i} = ['t = ',num2str(t(Nt_temp))];
    Nt_temp = Nt_temp + delta_Nt;
end
hold off
legend(name_cell);
xlabel('x');
ylabel('u'); 
title('Lax Scheme')

%% Lax-Wendroff Scheme
% \frac{u_{j}^{n+1}-u_{j}^{n}}{\Delta t}+\frac{u_{j+1}^{n}-u_{j-1}^{n}}{2\Delta x}=\frac{1}{2}\Delta t\frac{u_{j+1}^{n}-2u_{j}^{n}+u_{j-1}^{n}}{\Delta x^{2}}
% basic parameters
xspan = [0 3];
tspan = [0 3];
Nx = 100;
Nt = 85;
plots = 5;

% initial variables
delta_x = (xspan(2) - xspan(1))/Nx;
delta_t = (tspan(2) - tspan(1))/Nt;
c = delta_t / delta_x;

fprintf('Lax-Wendroff Scheme: c = %f\n',c);
x = zeros(1,Nx+1);
t = zeros(1,Nt+1);
u = zeros(Nt+1,Nx+1);

% initial mesh nodes
x = xspan(1):delta_x:xspan(2);   % x(i) spatial nodes
t = tspan(1):delta_t:tspan(2);   % t(i) time nodes

u(1,:) = sin(2*pi*(x));          % u(i,0) initial u using initial condition

% forward in time step
for n = 2:Nt+1
    for i = 2:Nx
        u(n,i) = (1-c^2)*u(n-1,i) + 0.5*c*(c-1)*u(n-1,i+1) + 0.5*c*(c+1)*u(n-1,i-1);  % inner points
    end
    u(n,1) = (1-c^2)*u(n-1,1) + 0.5*c*(c-1)*u(n-1,2) + 0.5*c*(c+1)*u(n-1,Nx);      % boundary points
    u(n,Nx+1) = (1-c^2)*u(n-1,Nx+1) + 0.5*c*(c-1)*u(n-1,2) + 0.5*c*(c+1)*u(n-1,Nx);      % boundary points
end
    
% plot results
% surf
f1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',f1);
colormap('cool');
[X,T] = meshgrid(x,t);
surf(T,X,u,'Parent',axes1,'EdgeAlpha',0.8,'EdgeColor','none');
box on;
grid on;
xlabel('t');
ylabel('x');
zlabel('u');
set(axes1,'FontSize',14,'LineWidth',2);
title('Lax-Wendroff Scheme');
% curves
figure('Color',[1 1 1]);
colormap jet

Nt_temp = 1;
name_cell = cell(1,plots);
for i = 1:plots
    delta_Nt = floor((Nt+1) / plots);
    plot(x,u(Nt_temp,:));
    hold on
    name_cell{i} = ['t = ',num2str(t(Nt_temp))];
    Nt_temp = Nt_temp + delta_Nt;
end
hold off
legend(name_cell);
xlabel('x');
ylabel('u'); 
title('Lax-Wendoff Scheme')

%% Leap-frog Scheme
% \frac{u_{j}^{n+1}-u_{j}^{n-1}}{2\Delta t}=\frac{u_{j+1}^{n}-u_{j-1}^{n}}{2\Delta x}
% basic parameters
xspan = [0 3];
tspan = [0 3];
Nx = 100;
Nt = 90;
plots = 5;

% initial variables
delta_x = (xspan(2) - xspan(1))/Nx;
delta_t = (tspan(2) - tspan(1))/Nt;
c = delta_t / delta_x;

fprintf('Leap-frog Scheme: c = %f\n',c);
x = zeros(1,Nx+1);
t = zeros(1,Nt+1);
u = zeros(Nt+1,Nx+1);

% initial mesh nodes
x = xspan(1):delta_x:xspan(2);   % x(i) spatial nodes
t = tspan(1):delta_t:tspan(2);   % t(i) time nodes

u(1,:) = sin(2*pi*(x));          % u(i,1) initial u using initial condition
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
    
% plot results
% surf
f1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',f1);
colormap('cool');
[X,T] = meshgrid(x,t);
surf(T,X,u,'Parent',axes1,'EdgeAlpha',0.8,'EdgeColor','none');
box on;
grid on;
xlabel('t');
ylabel('x');
zlabel('u');
set(axes1,'FontSize',14,'LineWidth',2);
title('Leap-frog Scheme');
% curves
figure('Color',[1 1 1]);
colormap jet

Nt_temp = 1;
name_cell = cell(1,plots);
for i = 1:plots
    delta_Nt = floor((Nt+1) / plots);
    plot(x,u(Nt_temp,:));
    hold on
    name_cell{i} = ['t = ',num2str(t(Nt_temp))];
    Nt_temp = Nt_temp + delta_Nt;
end
hold off
legend(name_cell);
xlabel('x');
ylabel('u'); 
title('Leap-frog Scheme');
%% end
