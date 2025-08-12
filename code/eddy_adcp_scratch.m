% scratch to create a mock data set of beam velocity given an eddy
% advecting by

% beam angle and width
phi = 25;
dphi = 3;

% cell boundaries
dist = 0.1:0.02:1.8;

% rankine vortex
R = 0.1; % vortex radius
v_max = 0.05; % m/s

%% test
x0 = 0;
y0 = 0;
t = 0;
U = 0.02;
V = 0;

x = linspace(-1,1,101);
y = linspace(-1,1,3);

vel = nan(length(x),length(y),2);
for i = 1:length(x)
    for j = 1:length(y)
        [uij,vij] = rankine(R,v_max,U,V,x0,y0,t,x(i),y(j));
        vel(i,j,1) = uij;
        vel(i,j,2) = vij;
    end
end

figure(99); clf; hold on
pcolor(x-0.5*diff(x(1:2)),y-0.5*diff(y(1:2)),vecnorm(vel,2,3)')
quiver(x,y,vel(:,:,1)',vel(:,:,2)','k')
axis tight


%%% subfunctions %%%
function [u,v] = rankine(R,v_max,U,V,x0,y0,t,x,y)
    C = R*v_max;
    % center of vortex
    xc = x0+U*t;
    yc = y0+V*t;

    phi = atan2(y-yc,x-xc);
    r = sqrt((y-yc)^2+(x-xc)^2);

    if r<=R
        v_phi = C*r/R^2;
    else
        v_phi = C/r;
    end

    u = -v_phi*sin(phi);
    v = v_phi*cos(phi);
end