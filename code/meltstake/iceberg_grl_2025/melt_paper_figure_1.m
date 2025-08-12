% Script to create figure 1 for the melt rate paper. This is a
% multi-panel figure that showcases the relevant data streams coming from
% the meltstake and how we calculate forcing. Redux of former figure
% versions.
%
% KJW
% 4 Sep 2024

clear

% choose deployment
dep_num = 8;

ms_tbl = loadMSInfo(dep_num);
dep_name = ms_tbl.Folder{1};

tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_scales = readtable(tbl_path,'sheet','manualwindows');
idx_tbl = ms_scales.Number == dep_num;
ms_scales = ms_scales(idx_tbl,:);
[T_scale,T_ci] = msTable2Vector(ms_scales.T);
[u_scale,u_ci] = msTable2Vector(ms_scales.u2);
[S_scale,S_ci] = msTable2Vector(ms_scales.S);
[m,m_ci] = msTable2Vector(ms_scales.m);

t_wnds = ms_scales.Start + minutes(0.5*ms_scales.Duration);

% load data
raw_path = fullfile('F:meltstake\data\raw',dep_name);
proc_path = fullfile('F:meltstake\data\proc',dep_name);
% T
load(fullfile(raw_path,'rbr','T.mat'))
% S
load(fullfile(proc_path,'ctd.mat'))
% adcp
load(fullfile(raw_path,'adcp','adcp.mat'))
% adv
load(fullfile(raw_path,'adv','adv.mat'),'adv')
load(fullfile(proc_path,'adv','pck_edges.mat'))

%% process and transform ADCP
cor_min = 60;
% adcp.burst.vel(adcp.burst.cor<cor_min) = nan;
% adcp.burst.vel = unwrapVel(adcp.burst.vel,.2,[],32);
adcp = msADCPTransform(adcp,cor_min);

% adcpQuickPlot(figure(1001),adcp,'vel_ice',.15*[-1 1],NaN,[0 1],1)

%% plottable data
%%% T %%%
T_mean = T(1).values;
for i = 2:3
    T_mean = T_mean + T(i).values;
end
T_mean = T_mean/3;
t_T = dn2dt(T(1).time);

%%% S %%%
ms_depth = str2double(ms_tbl.Depth_m_{1});
idxz_S = abs(ctd.z-ms_depth)<=1.5;
S = mean(ctd.S(idxz_S,:),'omitnan');

%%% vel %%%
vel = adcp.burst.vel_ice;
vel_scale = adcp.burst.vel_ice;
vel_smth = adcp.burst.vel_ice;

% smooth
% spatial
for k = 1:size(vel,1)
    for p = 1:size(vel,3)
        vel_smth(k,:,p) = hannFilter(squeeze(vel_smth(k,:,p)),3);
    end
end
vel_plot = vel_smth;

% calculate magnitude
vel_mag_smth = vecnorm(vel_smth,2,3);
% profile max and mean
vel_smth_max = max(vel_mag_smth,[],2);

% time
hann_dt = 1*30; % s
hann_width = ceil(hann_dt*adcp.burst.samplerate);

vel_smth_max = hannFilter(vel_smth_max,hann_width);

for i = 1:size(vel,2)
    for j = 1:size(vel,3)
        vel(:,i,j) = hannFilter(vel(:,i,j),hann_width);
        vel_scale(:,i,j) = hannFilter(vel_scale(:,i,j),8);
    end
end

% 
vel_vol = squeeze(mean(vel_scale(:,adcp.burst.range<=0.5,:),2,'omitnan')); % volume average of u,w,v

% downsample
idxt = ceil(hann_width/2):ceil((hann_width/10)):size(vel,1);
t_vel = dn2dt(adcp.burst.time(idxt));

idxr = adcp.burst.range <= 1;
r_vel = adcp.burst.range(idxr);
r_vel_ice = 0.65-r_vel;

vel_u = vel(idxt,idxr,1);
vel_w = vel(idxt,idxr,2);
vel_w(:,r_vel>0.65) = 0.5*vel_w(:,r_vel>0.65);

%%% adv ice %%%
% choose beam
bn = getBeamNumbers(ms_scales);

% downsample
rlim = [100 290];
idxr = adv.pck_dist(2,:)>=rlim(1) & adv.pck_dist(2,:)<=rlim(2);
r_ice = adv.pck_dist(2,idxr);
t_ice = adv.pck_time;
amp_ice = adv.pck_amp(:,idxr,bn);

% ice position
pos_ice = medianFilter(edges.pos(:,bn),3);

% fix drill jump
for i = 1:length(tdrill)
    drill_idx = find(edges.time==tdrill(i));
    drill_depth = pos_ice(drill_idx+2) - pos_ice(drill_idx);
    pos_ice(drill_idx+[0 1]) = nan;
    
    pos_ice2 = pos_ice;
    pos_ice2(edges.time<tdrill(i)) = pos_ice2(edges.time<tdrill(i)) + drill_depth - 2;
end

if isempty(tdrill)
    pos_ice2 = pos_ice;
end

%% plot
fs = 11;
lw = 1;
clear ax

% create and size figure
figsize = [10 16];
fig = figure(1); clf
setFigureSize(fig,figsize);

% padding
pxi = .015;
pxo = .2;
pyi = .01;
pyo = .05;
pad = [pxi pyi pxo pyo];
shift = [0 .01];

% panel labels
% panel_lbls = {'$u$','$w$','max($|\vec{u}|$)','$T_w$','$S_w$','ADV echogram'};
panel_lbls = {'u','w','','','','ADV echogram'};

% temperature things
T_clr = [0.2 1 0.5]/1.2;
T_clrs = T_clr'*[1 0.5 0];
T_lgd_lbls = {'near','mid','far'};

% vel things
vel_lim = 0.151*[-1 1];

% u
ax(1) = axes(fig,'position',axgridpos(7,1,1,pad,shift));
pcolor(t_vel,r_vel_ice,vel_u')
shading flat
cmocean('bal')
clim(vel_lim)
set(ax(1),'xticklabel',{})
ylabel('y [m]','fontsize',fs)
ylim([-0.09 max(r_vel_ice)])
ax(1).FontSize = fs-2;

% w
ax(2) = axes(fig,'position',axgridpos(7,1,2,pad,shift));
pcolor(t_vel,r_vel_ice,vel_w')
shading flat
cmocean('bal')
clim(vel_lim)
set(ax(2),'xticklabel',{})
ylabel('y [m]','fontsize',fs)
ylim([-0.09 max(r_vel_ice)])
ax(2).FontSize = fs-2;

cbar1 = colorbar;
cbar1.Position = cbarpos(ax(2),.005,.025);
cbar1.Position(4) = 2*cbar1.Position(4) + pyi;
cbar1.FontSize = fs-3;
cbar1.Label.String = 'm/s';
cbar1.Label.FontSize = fs;

% (u2+w2)^1/2
ax(3) = axes(fig,'position',axgridpos(7,1,3,pad,shift),'fontsize',fs-2);
hold on
plot(dn2dt(adcp.burst.time),vel_smth_max,'linewidth',1,'color',0.*[1 1 1])
ylabel('U_{max} [m/s]','fontsize',fs)
set(ax(3),'xticklabel',{})
box on

% T
ax(4) = axes(fig,'position',axgridpos(7,1,4,pad,shift),'fontsize',fs-2);
hold on
clear T_handles
for i = 1:length(T)
    T_handles(i) = plot(dn2dt(T(i).time),hannFilter(T(i).values,2*5),'color',T_clrs(:,i));
end
box on
ylabel('T_w [\circC]','fontsize',fs)
set(ax(4),'xticklabel',{})
axpos = ax(4).Position;
ax(4).YTick = 4:2:10;

% lgd = legend(ax(4),T_lgd_lbls,'orientation','vertical','location','southeastoutside',...
%     'fontsize',fs-2,'edgecolor','k');
% ax(4).Position = axpos;
% lgd.Position = lgd.Position - [.15 .03 0 0];

% S
ax(5) = axes(fig,'position',axgridpos(7,1,5,pad,shift),'fontsize',fs-2);
hold on
plot(ctd.t0,S,'o','linewidth',1,'color',0.1*[1 1 1],'markersize',3,'markerfacecolor',colors(2))
ylabel('S_w [psu]','fontsize',fs)
set(ax(5),'xticklabel',{})
box on
ylim([19. 28.5])
ax(5).YTick = 20:3:28;

lgd = legend(ax(5),T_handles,T_lgd_lbls,'orientation','vertical','location','southeastoutside',...
    'fontsize',fs-2,'edgecolor','k');
ax(4).Position = axpos;
lgd.Position = lgd.Position + [.15 .09 0 0];

% ice position
ax(6) = axes(fig,'position',axgridpos(7,1,7,pad,shift),'fontsize',fs-2);
ax(6).Position(4) = ax(5).Position(2) - ax(6).Position(2) - pyi;
hold on
pcolor(t_ice,r_ice/10,amp_ice')
shading flat
plot(edges.time,pos_ice/10,'k-','linewidth',lw)
plot(edges.time,pos_ice2/10,'k-')
for i = 1:length(t_wnds)
    tlim = [ms_scales.Start(i) ms_scales.End(i)];
    idxti = edges.time >= tlim(1) & edges.time <= tlim(2);
    ti = edges.time(idxti);
    ti = ti([1 end]) + minutes(2)*[1 -1]';

    pos_mean = mean(pos_ice2(idxti),'omitnan');
    plot(ti,0.1*pos_mean+m(i)*hours((ti-t_wnds(i)))-1,...
        'w--','linewidth',1,'markersize',8)
end
ylim(rlim/10)
clim([60 180])
ylabel({'distance','from ADV [cm]'},'fontsize',fs)
cbar2 = colorbar;
cbar2.Position = cbarpos(ax(6),.005,.025);
cbar2.Label.String = 'counts';
cbar2.Label.FontSize = fs;
% cbar2.Label.Interpreter = 'latex';

% window bars at top
ax(7) = axes(fig,'position',ax(1).Position);
hold on
for i = 1:length(t_wnds)
    plot([ms_scales.Start(i) ms_scales.End(i)]+minutes(1)*[1 -1],1.1*[1 1],'k|-',...
        'linewidth',1,'markersize',6)
    text(t_wnds(i),1.1,sprintf('segment %d',i),'interpreter','none','HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',fs-1)
end
ylim([0 1])
set(ax(7),'clipping','off')
set(ax(7),'visible','off')

for i = 1:6
%     text(ax(i),.01,.995,sprintf('(%c) %s',char(97+i),panel_lbls{i}),'units','normalized','fontsize',fs,...
%         'verticalalignment','top','horizontalalignment','left')
    text(ax(i),.01,.995,sprintf('%s',panel_lbls{i}),'units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left')
end

linkaxes(ax,'x')
xlim([ms_scales.Start(1) ms_scales.End(end)])

%% sideways!!!
% % create and size figure
% figsize = [17 8];
% fig = figure(2); clf
% setFigureSize(fig,figsize);
% clear ax
% 
% % padding
% pad = [.01 .01 .1 .15];
% shift = [.0 .0];
% 
% % panel labels
% panel_lbls = {'$u$','$w$','max($|\vec{u}|$)','$T$','$S$','ADV echogram'}; 
% 
% % temperature things
% T_clr = [0.2 1 0.5]/1.2;
% T_clrs = T_clr'*[1 0.5 0];
% T_lgd_lbls = {'near','mid','far'};
% 
% % vel things
% vel_lim = 0.151*[-1 1];
% 
% % u
% ax(1) = axes(fig,'position',axgridpos(1,7,1,pad,shift),'fontsize',fs-2);
% pcolor(r_vel_ice,t_vel,vel_u)
% shading flat
% cmocean('bal')
% clim(vel_lim)
% set(ax(1),'yticklabel',{})
% xlabel({'ice','distance (m)'},'fontsize',fs)
% xlim([-0.09 max(r_vel_ice)])
% 
% % w
% ax(2) = axes(fig,'position',axgridpos(1,7,2,pad,shift),'fontsize',fs-2);
% pcolor(r_vel_ice,t_vel,vel_w)
% shading flat
% cmocean('bal')
% clim(vel_lim)
% set(ax(2),'yticklabel',{})
% xlabel({'ice','distance (m)'},'fontsize',fs)
% xlim([-0.09 max(r_vel_ice)])
% 
% cbar1 = colorbar;
% cbar1.Position = cbarpos(ax(2),.005,.025);
% cbar1.Position(4) = 2*cbar1.Position(4) + pyi;
% cbar1.Label.String = 'm/s';
% cbar1.Label.FontSize = fs;
% cbar1.Label.Interpreter = 'latex';
% 
% % (u2+w2)^1/2
% ax(3) = axes(fig,'position',axgridpos(1,7,3,pad,shift),'fontsize',fs-2);
% hold on
% plot(vel_smth_max,dn2dt(adcp.burst.time),'linewidth',1,'color',0.3*[1 1 1])
% xlabel('m/s','fontsize',fs)
% set(ax(3),'yticklabel',{})
% box on
% 
% % T
% ax(4) = axes(fig,'position',axgridpos(1,7,4,pad,shift),'fontsize',fs-2);
% hold on
% clear T_handles
% for i = 1:length(T)
%     T_handles(i) = plot(hannFilter(T(i).values,2*5),dn2dt(T(i).time),'color',T_clrs(:,i));
% end
% box on
% xlabel('$^\circ$C','fontsize',fs)
% set(ax(4),'yticklabel',{})
% axpos = ax(4).Position;
% 
% % S
% ax(5) = axes(fig,'position',axgridpos(1,7,5,pad,shift),'fontsize',fs-2);
% hold on
% plot(S,ctd.t0,'o','linewidth',1,'color',0.1*[1 1 1],'markersize',3,'markerfacecolor',colors(2))
% % errorbar(t_wnds,S_scale,S_ci,'ko','linewidth',1,'markerfacecolor',colors(3),'markersize',5)
% xlabel('psu','fontsize',fs)
% set(ax(5),'yticklabel',{})
% box on
% xlim([20 27.5])
% 
% lgd = legend(ax(5),T_handles,T_lgd_lbls,'orientation','vertical','location','southeastoutside',...
%     'fontsize',fs-2,'edgecolor','k');
% ax(4).Position = axpos;
% lgd.Position = lgd.Position + [.15 .09 0 0];
% 
% % ice position
% ax(6) = axes(fig,'position',axgridpos(1,7,6,pad,shift),'fontsize',fs-2);
% ax(6).Position(3) = 2*ax(6).Position(3) + pad(1);
% hold on
% pcolor(r_ice/10,t_ice,amp_ice)
% shading flat
% plot(pos_ice/10,edges.time,'k-','linewidth',lw)
% plot(pos_ice2/10,edges.time,'k-')
% for i = 1:length(t_wnds)
%     tlim = [ms_scales.Start(i) ms_scales.End(i)];
%     idxti = edges.time >= tlim(1) & edges.time <= tlim(2);
%     ti = edges.time(idxti);
%     ti = ti([1 end]) + minutes(2)*[1 -1]';
% 
%     pos_mean = mean(pos_ice2(idxti),'omitnan');
%     plot(0.1*pos_mean+m(i)*hours((ti-t_wnds(i)))-1,ti,...
%         'w--','linewidth',1,'markersize',8)
% end
% xlim(rlim/10)
% clim([50 180])
% xlabel({'distance','from ADV (cm)'},'fontsize',fs)
% xlabel({'distance from ADV (cm)'},'fontsize',fs)
% set(ax(6),'yticklabel',{})
% cbar2 = colorbar;
% cbar2.Position = cbarpos(ax(6),.005,.025);
% cbar2.Label.String = 'counts';
% cbar2.Label.FontSize = fs;
% cbar2.Label.Interpreter = 'latex';
% 
% % window bars at top
% ax(7) = axes(fig,'position',ax(1).Position + [0 0 ax(1).Position(3)+pad(1) 0]);
% hold on
% for i = 1:length(t_wnds)
%     plot(-0.1*[1 1],[ms_scales.Start(i) ms_scales.End(i)]+minutes(1)*[1 -1],'k_-',...
%         'linewidth',1,'markersize',6)
%     text(-0.1,t_wnds(i),sprintf('Segment %d',i),'interpreter','latex','HorizontalAlignment','center',...
%         'VerticalAlignment','bottom','fontsize',fs,'Rotation',90)
% end
% xlim([0 1])
% xlabel('ice distance (m)','fontsize',fs)
% set(ax(7),'clipping','off')
% set(ax(7),'visible','off')
% 
% for i = 1:6
%     text(ax(i),.01,.98,sprintf('(%c) %s',char(97+i),panel_lbls{i}),'units','normalized','fontsize',fs,...
%         'verticalalignment','top','horizontalalignment','left')
% end
% 
% linkaxes(ax,'y')
% ylim([ms_scales.Start(1) ms_scales.End(end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bn = getBeamNumbers(ms_scales)
    bn_strs = strsplit(ms_scales.Beam{1},',');
    bn = str2double(bn_strs{1});
end


