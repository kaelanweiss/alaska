% Script to plot ADV backscatter profiles for melt rate methods
%
% KJW
% 27 Jun 2026
clear

seg_tbl = loadMSInfo(28,'segments');

folder = 'ms01_20240717_0050';
load(fullfile('F:/meltstake/data/raw',folder,'adv','adv.mat'))
load(fullfile('F:/meltstake/data/proc',folder,'adv','pck_edges.mat'))
load('F:/adv/mean_profile.mat')

%% create backscatter plot data
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

%% fix the silly business in the third segment
idxt = edges.time>=seg_tbl.Start(3) & edges.time<=seg_tbl.End(3);
idx_bad = idxt & edges.pos(:,2)>200;
edges.pos(idx_bad,2) = interp1(edges.time(~idx_bad),edges.pos(~idx_bad,2),edges.time(idx_bad));
% edges.pos(idxt,2) = hannFilter(edges.pos(idxt,2),3);

%% plot one beam
fig2 = figure(2); clf
setFigureSize(fig2,[9 5.5]);

lw = 2;
fs = 11;

pad = [.03 .05 .16 .13];
shift = [-0.03 0.07];
clear ax
ax2 = axes(fig2,'position',axgridpos(1,1,1,pad,shift));

% plot
% t_min = minutes(edges.time-datetime(2024,7,17,0,52,0));%seg_tbl.Start(1));
t_min = edges.time;
hold(ax2,'on')
pcolor(ax2,t_min,adv.pck_dist(1,:)/10,amp(:,:,2)')
shading(ax2,'flat')
plot(ax2,t_min,(medianFilter(edges.pos(:,2),1)-4)/10,'w-','linewidth',lw)
set(ax2,'fontsize',fs)

ylim(ax2,[120 200]/10)

cbar2 = colorbar(ax2,'position',cbarpos(ax2,.01,.02));
cbar2.Label.String = 'normalized amplitude';
cbar2.Label.FontSize = fs-1;


% xlabel(ax2,'time since deployment [min]','fontsize',fs)
ylabel(ax2,'ADV distance [cm]','fontsize',fs)

clim(ax2,[1.05 1.25])
cmocean('hal',ax2)

% xlim(ax2,[20 60])
xlim(ax2,datetime(2024,7,17,0,52,0)+minutes([20 60]))