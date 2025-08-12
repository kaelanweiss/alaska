% Script to plot Spider boi
% KJW
% 9 Sep 2022
clear;

raw_path = 'F:/Alaska2022/data/iceberg_surveys/raw/';
mat_path = 'F:/Alaska2022/data/iceberg_surveys/mat/';
berg = '20220824_singingflower';
depname = fullfile(raw_path,berg,'spider');

% load deployment time
tlim = parseDeploymentTimes(depname);

% adcp
if exist(fullfile(mat_path,berg,'spider','adcp.mat'),'file')
    load(fullfile(mat_path,berg,'spider','adcp.mat'))
else
    adcp = parseNortekDeployment(depname,{},tlim(1),tlim(2));
end
qa_corr = adcp.burst.cor>=25;

% rbr
if exist(fullfile(mat_path,berg,'spider','rbr.mat'),'file')
    load(fullfile(mat_path,berg,'spider','rbr.mat'))
else
    rbr = parseRBRDeployment(depname,tlim(1),tlim(2));
end

%% velocity
flds = {'vel','cor','amp'};
ttls = {'right','up','left','down'};
cmaps = {'balance','matter','amp'};
units = {'m/s','%','dB'};
%clims = {0.08*[-1 1],[20 100],[53 85]};
clims = {0.04*[-1 1],[50 100],[60 85]};
iRmax = 40;
dt = 1;
ax = [];
fs = 14;
for i = 1:length(flds)
    figure(i); clf;
%     subplot(adcp.burst.nbeams+1,1,1);
%     ax(end+1) = gca;
%     hold on;
%     for j = 1:length(rbr)
%         plot(rbr(j).time,rbr(j).values(:,1));
%     end
%     xlim(tlim);
%     datetick('x','mmmdd HH:MM','keeplimits');
%     ylabel('temp [C]');
%     grid;
    for j = 1:adcp.burst.nbeams
        temp = adcp.burst.(flds{i})(1:dt:end,1:iRmax,j);
        if strcmp(flds{i},'vel')
            temp(~qa_corr(:,1:iRmax,j)) = NaN;
        end
        %subplot(adcp.burst.nbeams+1,1,j+1);
        subplot(adcp.burst.nbeams,1,j);
        ax(end+1) = gca;
        pcolor(adcp.burst.time(1:dt:end),adcp.burst.range(1:iRmax),temp');
        shading flat;
        ylabel('range (m)','fontsize',fs);
        if j == 1
            axpos1 = get(gca,'position');
        elseif j == adcp.burst.nbeams
            axpos2 = get(gca,'position');
            cbar = colorbar;
            cbar.Position = [sum(axpos2([1 3]))+.01 axpos2(2) .015 sum(axpos1([2 4]))-axpos2(2)];
            cbar.Label.String = units{i};
            cbar.Label.FontSize = fs;
        end
        clim(clims{i});
        xlim(tlim);
        datetick('x','HH:MM:SS','keeplimits');
        title(ttls{j},'fontsize',fs);
        cmocean(cmaps{i})
        set(gca,'YTick',0.2:0.3:2)
    end
end
linkaxes(ax,'x');

%% depth
figure(i+1); clf;
plot(adcp.attitude.time,adcp.attitude.pressure-0.9,'k-');
xlim(tlim);
datetick('x','HH:MM:SS','keeplimits');
grid;
ax(end+1) = gca;
linkaxes(ax,'x');
ylabel('pressure (dbar)','fontsize',fs)
        
        
        