% workspace for trying to figure out the best way to filter Spider ADCP
% data
%
% Levels of filtering:
%   1) do nothing
%   2) low pass filter
%   3) interpolate bad correlation values
%   4) interpolate bad correlation, low pass filter

addpath('../../../code/wengrove/')
addpath('../../../code/nash/')

try
    adcp;
catch
    load F:\alaska2022\data\iceberg_surveys\mat\20220824_teaparty\spider\adcp.mat
end

% set max range
rmax = 1.1;
[~,ir] = min(abs(adcp.burst.range-rmax));

% set time window
t1 = min(adcp.burst.time);
t2 = max(adcp.burst.time);
% wind_num = 5;
% t1 = windows(wind_num,1);
% t2 = windows(wind_num,2);

[~,it1] = min(abs(adcp.burst.time-t1));
[~,it2] = min(abs(adcp.burst.time-t2));

vel = adcp.burst.vel(it1:it2,1:ir,:);
cor = adcp.burst.cor(it1:it2,1:ir,:);
r = adcp.burst.range(1:ir);
t = adcp.burst.time(it1:it2);
dt = diff(t(1:2))*86400;

% plotting
vlim = 0.08;
fs = 12;
%% do nothing
plotADCP(1,t,r,vel,vlim,fs,'original data');

%% low-pass 1
vel_lp1 = vel;
f0 = 1/15;
for i = 1:size(vel,2)
    for j = 1:size(vel,3)
        vel_lp1(:,i,j) = lowpass(vel(:,i,j),f0,1/dt);
    end
end

plotADCP(2,t,r,vel_lp1,vlim,fs,sprintf('low-pass f = %.2f s^{-1} (matlab)',f0));

%% moving average 2
vel_lp2 = vel;
k = (2*floor((2/dt)/2))+1;
for i = 1:size(vel,2)
    for j = 1:size(vel,3)
        vel_lp2(:,i,j) = meanFilter(vel(:,i,j),k);
    end
end

plotADCP(3,t,r,vel_lp2,vlim,fs,sprintf('low-pass f = %.2f s^{-1} (Meagan), no interpolation',f0));

%% QA only
min_cor = 60;
qa = cor<min_cor;
vel_qa = vel;
vel_qa(qa) = nan;

plotADCP(4,t,r,vel_qa,vlim,fs,sprintf('qa only (correlation >= %d)',min_cor));


%% QA and interpolate
vel_int1 = interpolateADCP(vel,qa,8*3,1);

plotADCP(5,t,r,vel_int1,vlim,fs,sprintf('interpolation only'));

%% QA, interpolate, low-pass
vel_int_lp = interpolateADCP(vel,qa,4,1);
f0 = 1e-1;
for i = 1:size(vel,2)
    for j = 1:size(vel,3)
        vel_int_lp(:,i,j) = lowpass(vel_int_lp(:,i,j),f0,1/dt);%highFreqFilt(vel_int_lp(:,i,j),dt,f0);
    end
end

plotADCP(6,t,r,vel_int_lp,vlim,fs,sprintf('interpolation and low-pass f = %.2f s^{-1}',f0));


%%% subfunctions %%%
function [fig,axlist,cbar] = plotADCP(fnum,t,r,vel,vlim,fs,lbl_text)
    beams = {'right','up','left','down'};
    fig = figure(fnum); clf
    for i = 1:4
        axlist(i) = axes(fig,'position',axgridpos(4,1,i,0.08,0.05));
        pcolor(t,r,vel(:,:,i)')
        shading flat
        cmocean('bal')
        clim(vlim*[-1 1])
        ylabel('range (m)','fontsize',fs)
        title(beams{i},'fontsize',fs)
        datetick('x','keeplimits')
    end
    linkaxes(axlist)

    % colorbar
    pos1 = get(axlist(1),'position');
    pos4 = get(axlist(4),'position');
    cbar = colorbar;
    cbar.Position = [pos4(1)+pos4(3)+0.005 pos4(2) 0.01 pos1(2)+pos1(4)-pos4(2)];
    cbar.Label.String = ['\bf\leftarrow \rminward' repmat(' ',[1 10]) 'm/s' repmat(' ',[1 10]) 'outward \bf\rightarrow'];
    cbar.Label.FontSize = fs;

    % label text
    text(axlist(1),min(xlim),max(ylim)+0.15*diff(ylim),lbl_text,'FontSize',fs)
end
