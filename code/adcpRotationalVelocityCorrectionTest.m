% looking at the effect of the center of rotation of the ROV differing from
% the radial center of the ADCP beams on observed along-beam velocities

clear;

% rotation rate
f = 10; % deg/s
w = 2*pi*f/360; % rad/s

% distance between beam-center and center of rotation
d = .2; % m

% distance from center of rotation
r = 0.2:.2:60; % m

% beam angle
theta = 25; % deg

% calculate (just a bunch of trig)
sindphi = d.*sind(theta)./sqrt(r.^2+(r-d).^2*tand(theta)^2);
dv = w*r.*sindphi;
dv_limit = w*d*sind(theta)*cosd(theta);

fprintf('Rotation rate: %.2f deg/s\n',f);
fprintf('Largest velocity correction: %.4f m/s\n',max(dv));

% plot
figure(201); clf; hold on;
plot(r,dv,'.-');
plot(r([1 end]),r([1 end])*0+dv_limit,'k--');
xlabel('r');
ylabel('dv');
grid;
ylim([0 1.01*max(dv)]);

figure(202); clf;
plot(r,(dv-dv_limit)/dv_limit,'.-');
xlabel('r')
ylabel('fractional error');
grid;


