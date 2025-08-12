% Script to analyze how much pitching due to melting will bias the melt
% rate high. Based on better scratch paper geometry.
%
% KJW
% 4 Apr 2025

% Definitions (lengths in m)
%   phi: pitch angle due to drills melting
%   L: vertical distance between ADV head and drill axis
%   D: distance between screw tip and ADV head, measured along drill axis
%   beta: ice angle, 0 is vertical, positive is undercut
%   b: horizontal distance between screw tip and ice face, measured along
%       drill axis
%   theta (x): phi - beta, pitch in ice-normal coordinates (makes math easier)

clear

L = .156*cosd(30);
D = .15+.17; % screw-to-rod + ADV-offset

pc_tbl = loadMSInfo(26:28,'pitch_correction');
ms_tbl = loadMSInfo(26:28,'manualwindows');

% distance (d) and pitch derivative of distance (dd_dtheta), take theta (x)
% and b as input
dvec_fxn = @(x,b) [(D-b).*cos(x) - L*sin(x); (D-b).*sin(x) - L*sin(x).*tan(x)];
d_fxn = @(x,b) sqrt((D-b).^2 + L^2*tan(x).^2 - 2*(D-b).*L.*tan(x));
dd_dx_fxn = @(x,b) L./cos(x).^2.*sign(L*tan(x) - (D-b));

b_fxn = @(x,d) D - L*tan(x) - d; % negative soln puts wall between ADV and drill tip

%% sanity check
xpad = 0;
x = linspace(-pi/4+xpad,pi/4-xpad,101);
b_test = .15;

d_test = d_fxn(x,b_test);
dd_dx_test = dd_dx_fxn(x,b_test);

dd_dx_num = gradient(d_test)/diff(x(1:2));

figure(1); clf
plot(x,d_test)
xlabel('\theta')
ylabel('dist [m]')

figure(2); clf
hold on
plot(x/pi*180,dd_dx_test,'linewidth',1)
plot(x/pi*180,dd_dx_num,'rx')
xlabel('\theta')
ylabel('dd/d\theta [m]')
%% calculate melt bias for each segment
% measured quantities
d_adv = msTable2Vector(pc_tbl.adv_dist)/100;
[p_adv,dp_adv] = msTable2Vector(pc_tbl.pitch_mean);
[beta_adv,dbeta_adv] = msTable2Vector(pc_tbl.ice_slope);
th_adv = p_adv - beta_adv;
[pr_adv,dpr_adv] = msTable2Vector(pc_tbl.pitch_rate); % deg/hr

% calculated quantities
b_adv = b_fxn(th_adv*pi/180,d_adv);
d_calc = d_fxn(th_adv*pi/180,0.15);
b0 = 0.15;
dd_dx = dd_dx_fxn(th_adv*pi/180,b0);
m_bias = dd_dx.*pr_adv*(pi/180)*24; % convert deg to rad, 1/hr to 1/day

% uncertainty
dm_fxn = @(p,beta,pr,dp,dbeta,dpr) sqrt(2*(2*L*pr.*sin(p-beta)./cos(p-beta).^3).^2.*(dp.^2 + dbeta.^2) + (L^2./cos(p-beta).^4).*dpr.^2);
dm = dm_fxn(p_adv*pi/180,beta_adv*pi/180,pr_adv*(pi/180)*24,dp_adv*pi/180,dbeta_adv*pi/180,dpr_adv*(pi/180)*24);

% print (cm/hr)
for i = 1:length(m_bias)
    fprintf('%.2f (%.2f)\n',round(m_bias(i)*100/24,2),round(dm(i)*100/24,2))
end