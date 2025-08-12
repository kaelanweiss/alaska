% Script to create figure 1 for the melt rate thesis. This is a
% multi-panel figure that showcases the relevant data streams coming from
% the meltstake and how we calculate forcing.
%
% KJW
% 30 Jan 2024
set(0,'defaulttextinterpreter','latex');

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
cor_min = 80;
% adcp.burst.vel(adcp.burst.cor<cor_min) = nan;
% adcp.burst.vel = unwrapVel(adcp.burst.vel,.2,[],32);
adcp = msADCPTransform(adcp,cor_min);

adcpQuickPlot(figure(1001),adcp,'vel_ice',.1*[-1 1],NaN,[0 1],1)

%% plottable data
%%% T %%%
T_mean = T(1).values;
for i = 2:3
    T_mean = T_mean + T(i).values;
end
T_mean = T_mean/3;
t_T = dn2dt(T(1).time);

%%% S %%%
ms_depth = ms_tbl.Depth_m_;
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

%%
fs = 11;
lw = 1;
clear ax

% create and size figure
figsize = [12 12];
fig = figure(1); clf
setFigureSize(fig,figsize);
% set(fig,'units','centimeters'); 
% fpos = get(fig,'position');
% fpos(3:4) = figsize;
% set(fig,'position',fpos);
% set(fig,'paperunits','centimeters')
% set(fig,'papersize',figsize)

% padding
pxi = .015;
pxo = .25;
pyi = .015;
pyo = .1;

% panel labels
panel_lbls = {'max($|\vec{u}|$)','$T$','$S$','ADV echogram'}; 

% temperature things
T_clr = [0.2 1 0.5]/1.2;
T_clrs = T_clr'*[1 0.5 0];
T_lgd_lbls = {'near','mid','far'};

% (u2+w2)^1/2
ax(1) = axes(fig,'position',axgridpos(5,1,1,pxi,pyi,pxo,pyo));
hold on
plot(dn2dt(adcp.burst.time),vel_smth_max,'linewidth',1,'color',0.3*[1 1 1])
ylabel('m/s','fontsize',fs)
set(ax(1),'xticklabel',{})
box on

% T
ax(2) = axes(fig,'position',axgridpos(5,1,2,pxi,pyi,pxo,pyo));
hold on;
for i = 1:length(T)
    plot(dn2dt(T(i).time),hannFilter(T(i).values,2*5),'color',T_clrs(:,i))
end
box on
ylabel('$^\circ$C','fontsize',fs)
set(ax(2),'xticklabel',{})

axpos = ax(2).Position;
lgd = legend(ax(2),T_lgd_lbls,'orientation','vertical','location','southeastoutside',...
    'fontsize',fs-2,'edgecolor','k');
ax(2).Position = axpos;
lgd.Position = lgd.Position - [.01 -.008 0 0];

% S
ax(3) = axes(fig,'position',axgridpos(5,1,3,pxi,pyi,pxo,pyo));
hold on
plot(ctd.t0,S,'o','linewidth',1,'color',0.1*[1 1 1],'markersize',3,'markerfacecolor',colors(2))
% errorbar(t_wnds,S_scale,S_ci,'ko','linewidth',1,'markerfacecolor',colors(3),'markersize',5)
ylabel('psu','fontsize',fs)
set(ax(3),'xticklabel',{})
box on
ylim([20 27.5])

% ice position
ax(4) = axes(fig,'position',axgridpos(5,1,5,pxi,pyi,pxo,pyo));
ax(4).Position(4) = ax(3).Position(2) - ax(4).Position(2) - pyi;
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
%clim([50 180])
ylabel({'distance','from ADV (cm)'},'fontsize',fs)
cbar = colorbar;
cbar.Position = cbarpos(ax(4),.005,.025);
cbar.Label.String = 'counts';
cbar.Label.FontSize = fs;
cbar.Label.Interpreter = 'latex';

% window bars at top
ax(5) = axes(fig,'position',ax(1).Position);
hold on
for i = 1:length(t_wnds)
    plot([ms_scales.Start(i) ms_scales.End(i)]+minutes(1)*[1 -1],1.1*[1 1],'k|-',...
        'linewidth',1,'markersize',6)
    text(t_wnds(i),1.1,sprintf('Segment %d',i),'interpreter','latex','HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',fs)
end
ylim([0 1])
set(ax(5),'clipping','off')
set(ax(5),'visible','off')

for i = 1:4
    text(ax(i),.01,.98,sprintf('(%c) %s',char(98+i),panel_lbls{i}),'units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left')
end

linkaxes(ax,'x')
xlim([ms_scales.Start(1) ms_scales.End(end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bn = getBeamNumbers(ms_scales)
    bn_strs = strsplit(ms_scales.Beam{1},',');
    bn = str2double(bn_strs{1});
end


