% playing around with plane geometry, trying to find the plane defined by
% three points
%
% ax + by + z = z0
% or
% (ax + by + z)/z0 = 1

% get access to advZ2B
addpath ../../../../adv/code

% set up points
z_adv = [250-5 250+10 250];

A = nan(3);
for i = 1:3
    [~,adv_coords] = advZ2B(z_adv(i),i);
    A(i,:) = adv_coords;
end

% A = randi([-10 10],3); % rows are (x,y,z) points of plane
b = ones(3,1); %

% solve for unknown vector
x = A\b; % = (a b 1)'/z0

% solve for coefficients
z0 = 1/x(3);
a = x(1)*z0;
b = x(2)*z0;

% build normal vector
n = [a b 1]';
n = 50*n/norm(n);

% plot
% define ice plane
grid_limit = 1.2*max(abs(extrema(A)));
npts = 2;
[Xgrid,Ygrid] = meshgrid(grid_limit*linspace(-1,1,npts),grid_limit*linspace(-1,1,npts));
Z = z0 - a*Xgrid - b*Ygrid;

% define ADV receiver locations
rcvrs = nan(3);
for i = 1:3
    [~,adv_coords] = advZ2B(38,i);
    rcvrs(i,:) = adv_coords;
end

figure(1); clf; hold on
surf(Xgrid,Ygrid,Z)
for i = 1:3
    plot3([rcvrs(i,1) A(i,1)],[rcvrs(i,2) A(i,2)],[rcvrs(i,3) A(i,3)],'b')
end
plot3(A(:,1),A(:,2),A(:,3),'ro',markerfacecolor='r')
plot3([0 0],[0 0],extrema(Z),'k--')
plot3([0 n(1)],[0 n(2)],[z0 z0+n(3)],'k.-')
shading interp
xlabel('x')
ylabel('y')
zlabel('z')