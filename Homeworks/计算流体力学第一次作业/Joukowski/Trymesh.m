%% Mesh JUKOVSKI Airfoil
c = 0.5;
lambda = 0.1;
r = (1+lambda)*c:0.1:5;
theta = 0:pi/40:2*pi;
[R,THETA] = meshgrid(r,theta);
x = R.*cos(THETA) -lambda*c;
y = R.*sin(THETA);

W = x + sqrt(-1)*y;

% mapping
Z = 0.5*(W + c^2*(W.^(-1)) );
X = real(Z);
Y = imag(Z);

% plot
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
figure('Color',[1 1 1]);
% plot phi
for j = 1:jmax 
   for i = 1:imax-1
        plot([x(i,j),x(i+1,j)],[y(i,j),y(i+1,j)],'r.-')
        hold on;
   end
end

% plot varphi
for i = 1:imax 
   for j = 1:jmax - 1
        plot([x(i,j),x(i,j+1)],[y(i,j),y(i,j+1)],'b.-')
   end
end
axis equal

%% H-Zone 
H = 10;
h = 6 ;

% mesh by Cartesian
r = logspace(log10(10),log10(1000000),30)/10000;
theta = 0+eps:pi/20:pi-eps;
[R,THETA] = meshgrid(r,theta);
x = R.*cos(THETA);
y = R.*sin(THETA);

W = x + sqrt(-1)*y;

% mapping
b = H/h;
t = sqrt( (W-b^2) ./ (W-1) );
t1 = (1+t)./(1-t);
t2 = (b+t)./(b-t);

Z = H/pi*log(t1) - h/pi*log(t2);
X = real(Z);
Y = imag(Z);

% plot
[imax jmax] = size(X);
figure('Color',[1 1 1]);
% plot mesh in z plane
for j = 1:jmax 
   for i = 1:imax-1
        plot([X(i,j),X(i+1,j)],[Y(i,j),Y(i+1,j)],'r.-')
        hold on;
   end
end

% plot mesh in z plane
for i = 1:imax 
   for j = 1:jmax - 1
        plot([X(i,j),X(i,j+1)],[Y(i,j),Y(i,j+1)],'b.-')
   end
end

axis equal

figure('Color',[1 1 1]);
% plot phi
for j = 1:jmax 
   for i = 1:imax-1
        plot([x(i,j),x(i+1,j)],[y(i,j),y(i+1,j)],'r.-')
        hold on;
   end
end

% plot varphi
for i = 1:imax 
   for j = 1:jmax - 1
        plot([x(i,j),x(i,j+1)],[y(i,j),y(i,j+1)],'b.-')
   end
end
axis equal

pause
close all






