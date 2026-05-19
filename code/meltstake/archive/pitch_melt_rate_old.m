% Script to analyze how much pitching due to melting will bias the melt
% rate high. Based on some scratch paper geometry. This may end up in
% supporting info of the melt rate paper.
%
% KJW
% 22 Jan 2024

%clear

% Definitions (lengths in m)
%   phi: pitch angle due to drills melting
%   L: vertical distance between ADV head and drill axis
%   D: distance between screw tip and ADV head, measured along drill axis
%   beta: ice angle, 0 is vertical, positive is undercut
%   b: horizontal distance between screw tip and ice face

clear

n = 250;
phi_lim = [-10 10];
phi = linspace(phi_lim(1)*pi/180,phi_lim(2)*pi/180,n);
cp = cos(phi);
sp = sin(phi);
tp = tan(phi);

L = .156*cosd(30);
D = .2+.17;
beta = 0;
b = .15;

% Solve location of ADV head
r1 = [D*cp + L*sp;...
      D*sp - L*cp];

% Solve for intersection of ADV axis and ice face
if beta==0
    r2 = [b*ones(1,n);...
          (b - D*cp - L*sp).*tp + D*sp - L*cp];
else
    m1 = tp;
    m2 = tand(90 - beta);
    b1 = -tp.*(D*cp + L*sp) + D*sp - L*cp;
    b2 = -b*tand(90 - beta);
    
    r2 = [(b1 - b2)./(m2 - m1);...
          (b1.*m2 - b2.*m1)./(m2 - m1)];
end

% Calculate distance between ADV and ice
try
    d = vecnorm(r2-r1);
catch
    d = nan(1,n);
    for i = 1:n
        d(i) = norm(r2(:,i)-r1(:,i),2);
    end
end

%% Calculate the time derivative of d
% dd/dt = dd/dphi * dphi/dt
dphi_dt = -3.58*pi/(180*3600); % rad/s
dphi = diff(phi(1:2));
dd_dphi = [nan conv(d,[1 0 -1],'valid')/(2*dphi) nan];
dd_dt = dd_dphi.*dphi_dt*24*3600; % m/day

fprintf('%.3f m/day\n',max(abs(dd_dt)));

% distance from ice
dq = linspace(min(d),max(d),8);
phiq = interp1(d,phi,dq);
mq = interp1(phi,dd_dt,phiq);

% plot
figure(1); clf
hold on
plot(phi*180/pi,dd_dt)
xlabel('pitch (deg)')
ylabel('melt bias (m/day)')
title(sprintf('dphi/dt=%.2f deg/hr',dphi_dt*180/pi*3600))
plot(phiq*180/pi,mq,'k.')
for i = 2:length(dq)-1
    text(phiq(i)*180/pi,mq(i),sprintf('%.1f cm',round(100*dq(i),1)),'verticalalignment','bottom')
end
