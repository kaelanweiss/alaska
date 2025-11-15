% Script to compare velocity spectra between meltstakes and available
% mooring/other velocity measurements.
%
% KJW
% 26 Sep 2025
clear

addpath('..')
load glacier_clrs.mat

% collect velocity data
% meltstake
fprintf('Meltstake data...\n')
raw_dir = 'F:/meltstake/data/raw';
ms_tbl = loadMSInfo(26:28,'segments');
[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);

ms24(ndeps) = struct('dep_name',[],'time',[],'range',[],'vel_ice',[]);
for i = 1:ndeps
    % loading
    load(fullfile(raw_dir,dep_names{i},'adcp','adcp.mat'))
    % indexing
    idxd = strcmp(ms_tbl.Folder,dep_names{i});
    idxt = adcp.burst.time >= ms_tbl.Start(find(idxd,1)) & adcp.burst.time <= ms_tbl.End(find(idxd,1,'last'));
    idxr = adcp.burst.range <= min(ms_tbl.rmax(idxd));
    % grabbing data
    ms24(i).dep_name = dep_names{i};
    ms24(i).time = adcp.burst.time(idxt);
    ms24(i).range = adcp.burst.range(idxr);
    ms24(i).vel_ice = adcp.burst.vel_ice(idxt,idxr,[1 3 2]);
end
clear adcp

% 2024
% workhorse (N1)
fprintf('Workhorse 2024...\n')
wh24_file = 'Y:\proc\mooring\workhorse\adcp\SN14158_54.nc';
wh24 = loadMooringADCP(wh24_file);
% big (D1)
fprintf('Big 2024...\n')
% b24_1_file = 'Y:\proc\mooring\big\adcp\SN22479_01.nc';
% b24_1 = loadMooringADCP(b24_1_file);
b24_2_file = 'Y:\proc\mooring\big\adcp\SN22479_02.nc';
b24 = loadMooringADCP(b24_2_file);

% 2025
% workhorse
fprintf('Workhorse 2025...\n')
wh25_file = 'Z:\proc\mooring\workhorse\SN15645.nc';
wh25 = loadMooringADCP(wh25_file);

%% average/clean data
% meltstakes
for i = 1:ndeps
    ms24(i).uvw = squeeze(mean(ms24(i).vel_ice,2,'omitnan'));
    ms24(i).uvw = fixUVWNaN(ms24(i));
end

% wh24
idxt = wh24.time>=ms24(2).time(1) & wh24.time<=ms24(2).time(end);
idxr = wh24.distance>=22 & wh24.distance<=28;
wh24.uvw = squeeze(mean(cat(3,wh24.u(idxt,idxr),wh24.v(idxt,idxr),wh24.w(idxt,idxr)),2,'omitnan'));
wh24.uvw = fixUVWNaN(wh24);

% big24
idxt = b24.time>=(ms24(1).time(1)-hours(3)) & b24.time<=(ms24(1).time(end)+hours(6));
idxr = b24.distance>=30 & b24.distance<=80;
b24.uvw = squeeze(mean(cat(3,b24.u(idxt,idxr),b24.v(idxt,idxr),b24.w(idxt,idxr)),2,'omitnan'));
b24.uvw = fixUVWNaN(b24);

% wh25
idxt = wh25.time>=datetime(2025,7,9,0,0,0) & wh25.time<=datetime(2025,7,12,0,0,0);
idxr = wh25.distance>=10 & wh25.distance<=30; 
wh25.uvw = squeeze(mean(cat(3,wh25.u(idxt,idxr),wh25.v(idxt,idxr),wh25.w(idxt,idxr)),2,'omitnan'));
wh25.uvw = fixUVWNaN(wh25);

%% PSDs
% binning setup
nf = [4 8 16];
nf_stops = [0 1e-2 1e-1];

% melstakes
for i = 1:length(ms24)
    [ms24(i).P,ms24(i).f,ms24(i).Pb,ms24(i).fb,ms24(i).M,ms24(i).S_ci] = uvwPSD(ms24(i),nf,nf_stops);
end

% wh/big
[wh24.P,wh24.f,wh24.Pb,wh24.fb,wh24.M,wh24.S_ci] = uvwPSD(wh24,1*nf,nf_stops);
[b24.P,b24.f,b24.Pb,b24.fb,b24.M,b24.S_ci] = uvwPSD(b24,1*nf(1:2),nf_stops(1:2));
[wh25.P,wh25.f,wh25.Pb,wh25.fb,wh25.M,wh25.S_ci] = uvwPSD(wh25,2*nf,nf_stops);

%% figures
% big figure with all 6 data sources
clear ax
fig1 = figure(1); clf
pad = [.05 .08 .08 .08];
shift = [0.03 0];
for i = 1:6
    ax(i) = axes(fig1,'position',axgridpos(2,3,i,pad,shift));
    hold(ax(i),'on'); box(ax(i),'on'); grid(ax(i),'on')
    set(ax(i),'xtick',10.^[-6:2:2])
end
mr_clrs = [colors(1);colors(5);colors(7)];
% ms
h_ms = plotPSD(ax(1),ms24(1),vel_clrs);
plotPSD(ax(2),ms24(2),vel_clrs);
plotPSD(ax(3),ms24(3),vel_clrs);
% dep1
plotPSD(ax(4),ms24(1),0*vel_clrs+0.7);
h_mr = plotPSD(ax(4),b24,mr_clrs);
% dep2
plotPSD(ax(5),ms24(2),0*vel_clrs+0.7);
plotPSD(ax(5),wh24,mr_clrs);
% 2025
plotPSD(ax(6),ms24(3),0*vel_clrs+0.7);
plotPSD(ax(6),wh25,mr_clrs);

% axes
for i = 1:length(ax)
    xlim(ax(i),10.^[-5 1]);
    ylim(ax(i),10.^[-6 1]);
end

% legends
legend(ax(1),h_ms,{'u_{term}','v_{term}','w'},'fontsize',10,'location','northeast')
legend(ax(4),h_mr,{'u','v','w'},'fontsize',10,'location','northeast')

% titles
for i = 1:3
    t1 = ms24(i).time(1);
    t2 = ms24(i).time(end);
    title(ax(i),{sprintf('MS dep %d (%d m)',i,round(mean(ms_tbl.depth(ms_tbl.Number==25+i)))),sprintf('(%d-%d-%d %02d:%02d-%02d:%02d)',t1.Year,t1.Month,t1.Day,t1.Hour,t1.Minute,t2.Hour,t2.Minute)},'fontsize',12)
end
title(ax(4),'Big (D1) 2024','fontsize',12)
title(ax(5),'Workhorse (N1) 2024','fontsize',12)
title(ax(6),'Workhorse 2025','fontsize',12)

% labels
for i = [1 4]
    ylabel(ax(i),'Velocity PSD [m^2 s^{-2} cps^{-1}]','fontsize',12)
end
for i = 4:6
    xlabel(ax(i),'f [cps]','fontsize',12)
end

%% Functions
% load mooring data into structure
function data = loadMooringADCP(file)
    finfo = ncinfo(file);
    varnames = {finfo.Variables.Name}';
    data = struct();
    for i = 1:length(varnames)
        data.(varnames{i}) = ncread(file,varnames{i});
    end

    data.time = datetime(data.time,'convertfrom','epochtime');
end

% interpolate NaNs in mean adcp data
function uvw = fixUVWNaN(adcp)
    uvw = adcp.uvw;
    idx_nan = isnan(adcp.uvw);
    t = adcp.time;

    for i = 1:3
        ui_fix = interp1(t(~idx_nan(:,i)),uvw(~idx_nan(:,i),i),t(idx_nan(:,i)));
        uvw(idx_nan(:,i),i) = ui_fix;
    end

end

% calculate uvw psds
function [P,f,Pb,fb,M,S_ci] = uvwPSD(adcp,nf,nf_stops)
    % remove mean velocity
    uvw = adcp.uvw;
    for i = 1:3
        uvw(:,i) = uvw(:,i) - mean(uvw(:,i));
    end
    % calculate PSDs
    fs = seconds(mean(diff(adcp.time)))^-1;
    [Puu,f] = psd_matlab(uvw(:,1),fs,'notaper');
    f = f';
    P = repmat(Puu,[1,3]);
    for i = 2:3
        P(:,i) = psd_matlab(uvw(:,i),fs,'notaper');
    end
    % bin
    [Pb,fb,M,S_ci] = progressiveBin(f,P(:,1),nf,nf_stops);
    Pb = repmat(Pb,[1,3]);
    for i = 2:3
        Pb(:,i) = progressiveBin(f,P(:,i),nf,nf_stops);
    end
end

% plot PSD in axes
function h = plotPSD(ax,adcp,clrs)
    if isnan(clrs)
        clrs = [colors(1);colors(2);colors(3)];
    end
    for i = 3:-1:1
        h(i) = plot(ax,adcp.fb,adcp.Pb(:,i),'linewidth',1,'color',clrs(i,:));
    end
    set(ax,'xscale','log','yscale','log')
end

