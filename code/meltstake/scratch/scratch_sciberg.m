clear

load F:/AK_202307/adcp/proc/ms03_03_adcp.mat
load F:/AK_202307/adv/ms03/3/adv.mat
load F:/AK_202307/rbr/ms03/3/rbr.mat

%%
rmax = 0.5;
idxr = adcp.burst.range<=rmax;
vel = adcp.burst.vel_ice(:,idxr,1:3);
vel_rms = sqrt(squeeze(mean(vel.^2,2,'omitnan')));
vel_rms_all = sqrt(mean(vel_rms.^2,2,'omitnan'));
vel_scale = hannFilter(vel_rms_all,120*4+1);

%%
t_drill = datetime(2023,7,9,15,23,18);
dist = adv.dsvol_ave(1:end-1,:);
dist(dist==0) = nan;
dist(adv.burst_time>t_drill,:) = dist(adv.burst_time>t_drill,:)+80;
dist = max(dist,[],2);
% for i = 1:2
%     dist(:,i) = medianFilter(dist(:,i),1);
% end

dist = medianFilter(dist,5)/10;

t_adv = datenum(adv.burst_time);
plot(t_adv,dist)

%%
t_rbr = datenum(rbr(1).t);
T = hannFilter(mean([rbr(1).T rbr(2).T rbr(3).T],2),120*2+1);

%%
figure(100); clf
t = {adcp.burst.time, t_rbr, t_adv};
dat = {vel_scale, T, dist};
clr = {colors(1),colors(2),colors(4)};
YLIM = {[0 0.11],[4 12],[6 17.8]};
YLBL = {{'ocean velocity','m/s'},{'ocean temperature','^\circC'},{'ice distance (slope=melt rate)','cm'}};

clear ax
for i = 1:3
    ax(i) = axes(figure(100),'position',axgridpos(3,1,i,0.05,0.02,0.12,0.08));
    plot((t{i}-t{1}(1))*86400/60,dat{i},color=clr{i},linewidth=1)
    ylim(YLIM{i})
    ylabel(YLBL{i}{2},fontsize=12)
    text(ax(i),0.02,0.92,YLBL{i}{1},units='normalized',fontsize=12)
    if i<3
        set(ax(i),xticklabel={})
    else
        xlabel('time (min)',fontsize=12)
    end

end

hold(ax(3),'on')
plot(ax(3),[0 120],[7 11],'k--',linewidth=1)
plot(ax(3),[120 200],[10.5 17],'k--',linewidth=1)

