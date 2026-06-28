% Script to bin the boundary layer energetics by timescale to quantify how
% much the f>N bits contribute
%
% 27 Jun 2026
% KJW
clear

load ../glacier_grl_2025/platform_comparison_data/ms_psd.mat
load F:/meltstake/data/proc/glacier_2024_ctd.mat

for i = 1:length(adcp)
    % mean flow contribution
    vel_mean = squeeze(mean(mean(adcp(i).vel_ice,1,'omitnan'),2,'omitnan'));
    KE_mean = sum(vel_mean([1 3]).^2);

    % sort by N
    N = mean(ctd(i).N_ms);
    P_ice = sum(mean(adcp(i).P_ice(:,:,[1 3]),2),3);
    idx_slow = adcp(i).f <= N;
    idx_fast = adcp(i).f > N & adcp(i).f < 1;

    df = mean(diff(adcp(i).f));
    KE_slow = df*trapz(P_ice(idx_slow));
    KE_fast = df*trapz(P_ice(idx_fast));

    % fprintf('%d | %.5f | %.5f | %.5f | m^2/s^2\n',i,round(KE_mean,5),round(KE_slow,5),round(KE_fast,5))
    % fprintf('-----\n')
    % fprintf('%d & %.3f & %.3f & %.3f | m/s\n',i,round(sqrt(KE_mean),3),round(sqrt(KE_slow),3),round(sqrt(KE_fast),3))
    fprintf('%.3f & %.3f & %.3f\n',round(sqrt(KE_mean),3),round(sqrt(KE_slow),3),round(sqrt(KE_fast),3))
end
