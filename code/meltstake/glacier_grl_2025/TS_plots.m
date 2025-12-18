% Script to work through what info to include in the TS plots of the
% boundary layer characteristics figure. Also figure out how to get rid of
% unphysical TS trajectories from the sensor misbehaving.
%
% KJW
% 24 Nov 2025
clear

raw_dir = 'F:/meltstake/data/raw';
proc_dir = 'F:/meltstake/data/proc';

load platform_comparison_data\ms_psd.mat

dep_tbl = loadMSInfo(26:28);
seg_tbl = loadMSInfo(26:28,'segments');

%% load hobo, T, polly ctd data
clear hobo TS
S0 = [26.7 27.8 28];
T0 = [6 6.9 7.1];
ms_depth = [21 43 52]; % m
TS(3) = struct();
for i = 1:3
    % hobo
    hobo_load = load(fullfile(raw_dir,dep_tbl.Folder{i},'hobo','hobo.mat'));
    hobo(i) = hobo_load.hobo(2);

    % TS
    idxt = hobo(i).time >= dep_tbl.Start(i) & hobo(i).time <= dep_tbl.End(i);
    TS(i).t = hobo(i).time(idxt);
    TS(i).S = hannFilter(hobo(i).S(idxt),1);
    TS(i).T = hannFilter(hobo(i).T_cal(idxt),1);
%     TS(i).fm = (T0(3)-Ti)/(T0(3)+90); % melt fraction
%     TS(i).fs = 1-Ti/T0(i); % fresh fraction
%     TS(i).rho0 = gsw_rho(S0(i),T0(i),ms_depth(i)); % outer density
%     TS(i).rho = gsw_rho(Si,Ti,ms_depth(i)); % bl density
%     TS(i).buoy = (9.8/rho0)*(rho0-rhoi); 
end

%% adcp data
% where to sample for velocities
y_s = [0.21] + 1e-6*[-1 1]; % distance from wall 
for i = 1:length(adcp)
    adcp(i).range = round(double(adcp(i).range),4);
    y_max = adcp(i).range(end);
    idx_s = (y_max-adcp(i).range) >= min(y_s) & (y_max-adcp(i).range) <= max(y_s);
    adcp(i).idx_s = idx_s;
    adcp(i).uvw = squeeze(mean(adcp(i).vel_ice(:,idx_s,:),2));
end
if ~(sum(adcp(1).idx_s) == sum(adcp(2).idx_s) && sum(adcp(1).idx_s) == sum(adcp(3).idx_s))
    warning('Number of bins do not agree')
end

%% clean up unphysical TS trajectories
k_hann = 3;
d_phi = 10*pi/180; % rad
fig = figure(100); clf
for i = 1
    Si = hannFilter(TS(i).S,k_hann);
    Ti = hannFilter(TS(i).T,k_hann);
    plot(Si,Ti,'k.-')
    hold on
    dT = diff(hannFilter(TS(i).T,k_hann));
    dS = diff(hannFilter(TS(i).S,k_hann));
    phi = atan2(dT,dS);
    phi = phi + pi/2;
    phi(phi<0) = phi(phi<0)+2*pi;
    phi = phi - pi/2;
    idx = abs(phi)<d_phi | abs(phi-pi)<d_phi;
    [blocks,widths] = findBlocks(idx);
    plot(Si(idx),Ti(idx),'r.')

end

%% gridded TS
dS_grid = [.05, 1.5*.05, .05];
dT_grid = [.04, 1.5*.04, .04];

S_grid = cell(3,1);
T_grid = cell(3,1);
n_grid = cell(3,1);
umag_grid = cell(3,1);
w_grid = cell(3,1);

% overall structure
grid_TS(3) = struct('S',[],'T',[],'n',[],'u_mag',[],'w_prime',[]);

for i = 1:3
    % calculate grid
    S_grid{i} = (floor(min(TS(i).S/dS_grid(i))):ceil(max(TS(i).S/dS_grid(i))))'*dS_grid(i);
    T_grid{i} = (floor(min(TS(i).T/dT_grid(i))):ceil(max(TS(i).T/dT_grid(i))))'*dT_grid(i);

    % other fields
    u_mag = interp1(adcp(i).time,hannFilter(vecnorm(adcp(i).uvw,2,2),5),TS(i).t,'nearest');
%     w_prm = interp1(adcp(i).time,hannFilter(adcp(i).uvw(:,3),5),TS(i).t,'nearest');
    w_prm = interp1(adcp(i).time,adcp(i).uvw(:,3) - hannFilter(adcp(i).uvw(:,3),5*8*60),TS(i).t,'nearest');
    w_prm = interp1(adcp(i).time,adcp(i).uvw(:,3),TS(i).t,'nearest');
    w_prm = detrend(w_prm,'omitnan');

    % normalize
    w_prm = w_prm/std(w_prm,'omitnan');
    u_srt = sort(u_mag);
    u_90 = u_srt(floor(0.9*length(u_mag)));
    u_mag = u_mag/u_90;

    mean(w_prm,'omitnan')

    % grid counts
    nS = length(S_grid{i})-1;
    nT = length(T_grid{i})-1;
    for j = 1:nS
        S_lim = S_grid{i}(j:j+1);
        for k = 1:nT
            T_lim = T_grid{i}(k:k+1);
            % count
            idxjk = (TS(i).S>=S_lim(1) & TS(i).S<S_lim(2)) & (TS(i).T>=T_lim(1) & TS(i).T<T_lim(2));
            n_grid{i}(j,k) = sum(idxjk);
            % u magnitude
            umag_grid{i}(j,k) = mean(u_mag(idxjk),'omitnan');
            % w prime
            w_grid{i}(j,k) = mean(w_prm(idxjk),'omitnan');

            % get rid of underrepresented bins
            idx_clear = n_grid{i}<3;
            umag_grid{i}(idx_clear) = nan;
            w_grid{i}(idx_clear) = nan;
        end
    end

    % overall output
    grid_TS(i).S = S_grid{i};
    grid_TS(i).T = T_grid{i};
    grid_TS(i).n = n_grid{i};
    grid_TS(i).u_mag = umag_grid{i};
    grid_TS(i).w_prime = w_grid{i};
end

%% plot
fig = figure(10); clf
clear ax
pad = [.13 .08 .09 .1];
shift = [-.02 0];
for i = 1:3
    ax(3*(i-1)+1) = axes(fig,'position',axgridpos(3,3,3*(i-1)+1,pad,shift));
    pcolor(ax(3*(i-1)+1),S_grid{i}(1:end-1),T_grid{i}(1:end-1),log10(n_grid{i}'))
    shading flat
    cbar = colorbar(ax(3*(i-1)+1),'position',cbarpos(ax(3*(i-1)+1),.005,.015));
    cbar.Label.String = 'log_{10}(count)';
    cbar.Label.FontSize = 10;
    colormap(ax(3*(i-1)+1),'parula')

    ax(3*(i-1)+2) = axes(fig,'position',axgridpos(3,3,3*(i-1)+2,pad,shift));
    pcolor(ax(3*(i-1)+2),S_grid{i}(1:end-1),T_grid{i}(1:end-1),umag_grid{i}')
    shading flat
    cbar = colorbar(ax(3*(i-1)+2),'position',cbarpos(ax(3*(i-1)+2),.005,.015));
    cbar.Label.String = '|u|/u_{90%}';
    cbar.Label.FontSize = 10;
    cmocean('mat',ax(3*(i-1)+2))
    clim([0.3 1])

    ax(3*(i-1)+3) = axes(fig,'position',axgridpos(3,3,3*(i-1)+3,pad,shift));
    pcolor(ax(3*(i-1)+3),S_grid{i}(1:end-1),T_grid{i}(1:end-1),w_grid{i}')
    shading flat
    cbar = colorbar(ax(3*(i-1)+3),'position',cbarpos(ax(3*(i-1)+3),.005,.015));
    cbar.Label.String = 'w''/\sigma_w';
    cbar.Label.FontSize = 10;
    cmocean('curl',ax(3*(i-1)+3))
    clim(1.2*[-1 1])
end

for i = 1:length(ax)
    xlabel(ax(i),'S [psu]')
    ylabel(ax(i),'T [\circC]')
end
for i = 1:3
    ylabel(ax(3*(i-1)+1),{sprintf('deployment %d',i),'T [\circC]'})
end
title(ax(1),'measurement density')
title(ax(2),'velocity magnitude')
title(ax(3),'vertical velocity anomaly')
linkaxes(ax(1:3));linkaxes(ax(4:6));linkaxes(ax(7:9))
xlim(ax(1),[25 26.5]); ylim(ax(1),[4.5 5.8])
xlim(ax(4),[25 S_grid{2}(end-1)]); ylim(ax(4),[5.2 T_grid{2}(end-1)])
xlim(ax(7),[26 S_grid{3}(end-1)]); ylim(ax(7),[6 T_grid{3}(end-1)])

% save
% save platform_comparison_data\grid_TS.mat grid_TS

