% workspace for trying to figure out the best way to filter Spider ADCP
% data
%
% Levels of filtering:
%   1) do nothing
%   2) low pass filter

try
    adcp;
catch
    load F:\alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
end

% set max range
rmax = 0.9;
[~,ir] = min(abs(adcp.burst.range-rmax));

% set time window
t1 = min(adcp.burst.time);
t2 = max(adcp.burst.time);
% t1 = datenum(2022,8,24,21,50,0);
% t2 = datenum(2022,8,24,22,15,0);
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
vlim = 0.055;
fs = 12;
tlim = [datenum(2022,8,25,3,15,30) ...
        datenum(2022,8,25,3,42,30)];
%% do nothing
[~,ax1,~] = plotADCP(1,t,r,vel,tlim,vlim,fs,'original data');
% for i = 1:4
%     ax1(i).XTick = ax1(i).XTick(1:2:end);
% end
%% low-pass 1
vel_lp1 = vel;
f0 = 1/20;
for i = 1:size(vel,2)
    for j = 1:size(vel,3)
        vel_lp1(:,i,j) = lowpass(vel(:,i,j),f0,1/dt);
    end
end

[~,ax2,~] = plotADCP(2,t,r,vel_lp1,tlim,vlim,fs,sprintf('low-pass f = %.2f s^{-1}',f0));
% for i = 1:4
%     ax2(i).XTick = (floor(1440*ax2(i).XLim(1)):1:floor(1440*ax2(i).XLim(2)))/1440;
%     datetick(ax2(i),'x','mmmdd HH:MMZ','keeplimits','keepticks')
% end

%%% subfunctions %%%
function [fig,axlist,cbar] = plotADCP(fnum,t,r,vel,tlim,vlim,fs,lbl_text)
    beams = {'beam1: right','beam2: up','beam3: left','beam4: down'};
    fig = figure(fnum); clf
    for i = 1:4
        axlist(i) = axes(fig,'position',axgridpos(4,1,i,0.08,0.04));
        pcolor(t,r,vel(:,:,i)')
        shading flat
        cmocean('bal')
        clim(vlim*[-1 1])
        xlim(tlim)
        ylabel('range (m)','fontsize',fs)
        %title(beams{i},'fontsize',fs)
        text(1,1,beams{i},'fontsize',fs,'fontweight','bold','units','normalized',...
            'horizontalalignment','right','verticalalignment','top')
        datetick('x','mmmdd HH:MMZ','keeplimits')
    end
    linkaxes(axlist)

    % colorbar
    pos1 = get(axlist(1),'position');
    pos4 = get(axlist(4),'position');
    cbar = colorbar;
    cbar.Position = [pos4(1)+pos4(3)+0.005 pos4(2) 0.015 pos1(2)+pos1(4)-pos4(2)];
    cbar.Label.String = ['\bf\leftarrow \rminward' repmat(' ',[1 10]) 'm/s' repmat(' ',[1 10]) 'outward \bf\rightarrow'];
    cbar.Label.FontSize = fs;

    % label text
    text(axlist(1),0,1,lbl_text,'fontsize',fs,'units','normalized',...
        'verticalalignment','bottom')
end
