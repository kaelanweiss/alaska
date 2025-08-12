% Script to run spider adcp fields through an EOF routine
% 
% KJW
% 13 Feb 2023

%% load data
load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
load windows_singingflower_0824.mat

%% choose window, trim, qa
wnd = 1;
rmax = 0.8;
[~,irmax] = min(abs(adcp.burst.range-rmax));
idxt = adcp.burst.time>=windows(wnd,1) & adcp.burst.time<windows(wnd,2);

t = adcp.burst.time(idxt);
r = adcp.burst.range(1:irmax);
vel = adcp.burst.vel(idxt,1:irmax,:);
cor = adcp.burst.cor(idxt,1:irmax,:);

vel(cor<50) = nan;

% plot
figure(1); clf

