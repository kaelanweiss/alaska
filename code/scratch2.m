% looking at the effect of the center of rotation of the ROV differing from
% the radial center of the ADCP beams on observed along-beam velocities

clear;

% rotation rate
f = 12; % deg/s
w = 2*pi*f/360; % rad/s

% distance between beam-center and center of rotation
d = 0.1;%.0375; % m

% distance from center of rotation
r = 0.05:.005:2; % m

% beam angle
theta = 25; % deg

% calculate (just a bunch of trig)
sindphi = d.*sind(theta)./sqrt(r.^2+(r-d).^2*tand(theta)^2);
dv = w*r.*sindphi;

fprintf('Rotation rate: %.2f deg/s\n',f);
fprintf('Largest velocity correction: %.4f m/s\n',max(dv));

% plot
figure(1); clf; hold on;
plot(r,dv,'.-');

plot(r([1 end]),r([1 end])*0+w*d*sind(theta)*cosd(theta),'k--');

xlabel('r');
ylabel('dv');
grid;


