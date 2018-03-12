%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name: JoukowskiOMesh.m
%%%% Program Prupose: Generate Joukowski Airfoil O mesh
%%%% Aurthor : Yang Yang
%%%% Date : 2015.09.21
%%%% Version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Jukowski Airfoil parameters
c = 0.5;
lambda = 0.1;
%% Mesh parameters
N_r = 100;
N_theta = 200; 
r_span = [(1+lambda)*c,5];
theta_span = [0,2*pi]; 
%% Generate mesh for Jukowski Airfoil
r = meshfun(r_span,N_r);
theta = linspace(theta_span(1),theta_span(2),N_theta); 
[R,THETA] = meshgrid(r,theta);
x = R.*cos(THETA) -lambda*c;
y = R.*sin(THETA) ;
W = x + sqrt(-1)*y;

%% confroaml mapping
Z = 0.5*(W + c^2*(W.^(-1)) );
X = real(Z);
Y = imag(Z);

%% plot Mesh in Z plane
[imax jmax] = size(X);
figure('Color',[1 1 1]);
% plot phi
for j = 1:jmax 
   for i = 1:imax-1
        plot([X(i,j),X(i+1,j)],[Y(i,j),Y(i+1,j)],'r.-')
        hold on;
   end
end

% plot varphi
for i = 1:imax 
   for j = 1:jmax - 1
        plot([X(i,j),X(i,j+1)],[Y(i,j),Y(i,j+1)],'b.-')
   end
end

axis equal

%% Plot mesh in W plane
figure('Color',[1 1 1]);
% plot r
for j = 1:jmax 
   for i = 1:imax-1
        plot([x(i,j),x(i+1,j)],[y(i,j),y(i+1,j)],'r.-')
        hold on;
   end
end

% plot theta
for i = 1:imax 
   for j = 1:jmax - 1
        plot([x(i,j),x(i,j+1)],[y(i,j),y(i,j+1)],'b.-')
   end
end
axis equal

%% Output Mesh to tecplot file
fp = fopen('JoukowskiOMesh.dat','w');
fprintf(fp,'Title = JoukowskiOMesh\n');
fprintf(fp,'VARIABLES = "X", "Y"\n');
fprintf(fp,'ZONE I =%d, J =%d, F = point\n',imax,jmax);

for j = 1:jmax
   for i = 1:imax
        fprintf(fp,'%e, %e\n',X(i,j),Y(i,j));
   end
end
fclose(fp);