% Script to plot ADV backscatter profiles for melt rate methods
%
% KJW
% 11 Feb 2026
clear

seg_tbl = loadMSInfo(28,'segments');

folder = 'ms01_20240717_0050';
load(fullfile('F:/meltstake/data/raw',folder,'adv','adv.mat'))
load(fullfile('F:/meltstake/data/proc',folder,'adv','pck_edges.mat'))
load('F:/adv/mean_profile.mat')

%% create plot data
dt = seconds(diff(adv.pck_time));
idxt = true(length(adv.pck_time),1);

% trim unpaired scans
if dt(1) > dt(2)
    idxt(1) = false;
end
if dt(end-1) < dt(end)
    idxt(end) = false;
end
amp1 = adv.pck_amp(1:2:end,:,:);
amp2 = adv.pck_amp(2:2:end,:,:);
amp = (amp1+amp2)/2;
time = adv.pck_time(idxt);
time = time(1:2:end);

% normalize using background profile
amp = amp./repmat(mean_prof',[size(amp,1) 1 size(amp,3)]).^1;

% smooth
[nt,nc,nb] = size(amp);
for k = 1:nb
    ampk = amp(:,:,k);
    % smooth in space
    for j = 1:nt
        ampk(j,:) = hannFilter(ampk(j,:),5);
    end
    % smooth in time
    for j = 1:nc
        ampk(:,j) = hannFilter(ampk(:,j),1);
    end
    amp(:,:,k) = ampk;
end

%% filter edges by distance
% edges.pos(edges.pos(:,1)>180,1) = nan;
% edges.pos(edges.pos(:,2)>200,2) = nan;
% edges.pos(edges.pos(:,3)>175,3) = nan;

%% plot
fig = figure(1); clf
setFigureSize(fig,[16 6]);

lw = 1.5;
fs = 11;

pad = [.03 .05 .15 .12];
shift = [0 0.05];
clear ax
for i = 3:-1:1
    ax(i) = axes(fig,'position',axgridpos(1,3,i,pad,shift));
end

% plot
t_min = minutes(edges.time-datetime(2024,7,17,0,52,0));%seg_tbl.Start(1));
for i = 1:3
    hold(ax(i),'on')
    pcolor(ax(i),t_min,adv.pck_dist(1,:)/10,amp(:,:,i)')
    shading(ax(i),'flat')
    plot(ax(i),t_min,(medianFilter(edges.pos(:,i),3)-5)/10,'w-','linewidth',lw)
    % cmocean('ice',ax(i))
    % colormap(ax(i),'jet')
    set(ax(i),'fontsize',fs)
end

linkaxes(ax)
ylim(ax(1),[120 220]/10)

for i = 2:3
    set(ax(i),'yticklabels',{})
end

cbar = colorbar(ax(3),'position',cbarpos(ax(3),.01,.02));
cbar.Label.String = 'normalized amplitude';
cbar.Label.FontSize = fs;
cbar.Label.FontWeight = 'bold';

xlabel(ax(2),'time since deployment [min]','fontsize',fs,'fontweight','bold')
ylabel(ax(1),'distance [cm]','fontsize',fs,'fontweight','bold')

for i = 1:3
    clim(ax(i),[1.05 1.25])
end

% xlim(ax(1),[seg_tbl.Start(snum-1) seg_tbl.End(snum)])
% xlim(ax(1),[minutes(seg_tbl.Start(2)-seg_tbl.Start(1)) minutes(seg_tbl.End(3)-seg_tbl.Start(1))])
% xlim(ax(1),[])

%% plot one beam
fig2 = figure(2); clf
setFigureSize(fig2,[12 6]);

lw = 2;
fs = 12;

pad = [.03 .05 .15 .12];
shift = [-0.02 0.07];
clear ax
ax2 = axes(fig2,'position',axgridpos(1,1,1,pad,shift));

% plot
t_min = minutes(edges.time-datetime(2024,7,17,0,52,0));%seg_tbl.Start(1));
hold(ax2,'on')
pcolor(ax2,t_min,adv.pck_dist(1,:)/10,amp(:,:,2)')
shading(ax2,'flat')
plot(ax2,t_min,(medianFilter(edges.pos(:,2),3)-5)/10,'w-','linewidth',lw)
set(ax2,'fontsize',fs)

ylim(ax2,[120 200]/10)

cbar2 = colorbar(ax2,'position',cbarpos(ax2,.01,.02));
cbar2.Label.String = 'normalized amplitude';
cbar2.Label.FontSize = fs;
cbar2.Label.FontWeight = 'bold';

xlabel(ax2,'time since deployment [min]','fontsize',fs,'fontweight','bold')
ylabel(ax2,'distance [cm]','fontsize',fs,'fontweight','bold')

clim(ax2,[1.05 1.25])
cmocean('hal',ax2)

xlim(ax2,[20 60])