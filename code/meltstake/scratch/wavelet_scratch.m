clear

load ../glacier_grl_2025/platform_comparison_data/ms_psd.mat

[~,i0] = max(adcp(1).f'.*mean(adcp(1).P_ice(:,:,1),2));
f0 = adcp(1).f(i0);
T0 = 1/f0;

vel = squeeze(mean(adcp(1).vel_ice,2));

% band pass
df = mean(diff(adcp(1).f));
wpass = f0 + 1.5*df*[-1 1] - 10*df;
vel_f0 = bandpass(vel,wpass,8);