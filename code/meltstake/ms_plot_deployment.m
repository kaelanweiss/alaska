% Script to plot T, various u, and ice distance fields to find windows for
% calculating the melt parameterization. This is for manually prepping melt
% rates for OSM
%
% KJW
% 5 Feb 2024

clear

% choose deployment
dep_num = 14;

ms_tbl = loadMSInfo(dep_num);
dep_name = ms_tbl.Folder{1};

% load data
raw_path = fullfile('F:meltstake\data\raw',dep_name);
proc_path = fullfile('F:meltstake\data\proc',dep_name);
% T
load(fullfile(raw_path,'rbr','T.mat'))
% adcp
load(fullfile(raw_path,'adcp','adcp.mat'))
% adv
load(fullfile(raw_path,'adv','adv.mat'),'adv')
% adv ice edges
load(fullfile(proc_path,'adv','pck_edges.mat'))

% process and transform ADCP
adcp = msADCPTransform(adcp);

%% plottable data
%%% T %%%
T_mean = T(1).values;
T_all = repmat(T(1).values,[1 3]);
for i = 2:3
    T_mean = T_mean + T(i).values;
    T_all(:,i) = T(i).values;
end
T_mean = T_mean/2;
t_T = dn2dt(T(1).time);

%%% vel %%%
vel = adcp.burst.vel_ice;

% smooth
hann_dt = 5; % s
hann_width = ceil(hann_dt*adcp.burst.samplerate);

for i = 1:size(vel,2)
    for j = 1:size(vel,3)
        vel(:,i,j) = hannFilter(vel(:,i,j),hann_width);
    end
end

% downsample
idxt = ceil(hann_width/2):ceil((hann_width)):size(vel,1);
t_vel = dn2dt(adcp.burst.time(idxt));

idxr = adcp.burst.range <= 0.9;
r_vel = adcp.burst.range(idxr);

vel_u = vel(idxt,idxr,1);
vel_w = vel(idxt,idxr,2);

%%% adv ice %%%
% choose beam
bn = 1;

% downsample
rlim = [0 350];
idxr = adv.pck_dist(2,:)>=rlim(1) & adv.pck_dist(2,:)<=rlim(2);
r_ice = adv.pck_dist(2,idxr);
t_ice = adv.pck_time;
amp_ice = adv.pck_amp(:,idxr,bn);

% ice position
pos_ice = medianFilter(edges.pos(:,bn),3);
% drill_idx = find(edges.time==tdrill);
% drill_depth = pos_ice(drill_idx+2) - pos_ice(drill_idx);
% pos_ice(drill_idx+[0 1]) = nan;

% fix drill jump
% pos_ice2 = pos_ice;
% pos_ice2(edges.time<tdrill) = pos_ice2(edges.time<tdrill) + drill_depth - 2;



% plot
set(0,'defaulttextinterpreter','latex');
fs = 11;
lw = 1;

figure(42); clf
clear ax

pdox = 0.15;
pdoy = 0.1;
pdix = 0.05;
pdiy = 0.08;
shift = [-(pdox-0.05) 0 0 0];

% panel 1: Temperature
ax(1) = axes(figure(42),'position',axgridpos(3,2,1,pdix,pdiy,pdox,pdoy,'flip') + shift);
plot(t_T,T_all)
% hold on;
% for i = 1:length(T)
%     plot(dn2dt(T(i).time),T(i).values)
% end


% panel 2: u
ax(2) = axes(figure(42),'position',axgridpos(3,2,2,pdix,pdiy,pdox,pdoy,'flip') + shift);
pcolor(t_vel,r_vel,vel_u')
shading flat
clim(0.08*[-1 1])
cmocean('bal')

% panel 3: w
ax(3) = axes(figure(42),'position',axgridpos(3,2,3,pdix,pdiy,pdox,pdoy,'flip') + shift);
pcolor(t_vel,r_vel,vel_w')
shading flat
clim(0.08*[-1 1])
cmocean('bal')


% panel 4: ice distance
ax(4) = axes(figure(42),'position',axgridpos(3,2,4,pdix,pdiy,pdox,pdoy,'flip') + shift);
hold on
pcolor(t_ice,r_ice,amp_ice')
shading flat
plot(edges.time,pos_ice,'k-','linewidth',lw)
%plot(edges.time,pos_ice2,'k-')
ylim(rlim)
%cmocean('rai')
clim([50 180])

% panel 5: 
ax(5) = axes(figure(42),'position',axgridpos(3,2,5,pdix,pdiy,pdox,pdoy,'flip') + shift);
plot(t_T,T(1).values);

% panel 6: 
pos6 = axgridpos(3,2,6,pdix,pdiy,pdox,pdoy,'flip') + shift;
ax(6) = axes(figure(42),'position',pos6);
plot(t_T,T(1).values);

% panel 7: 
expand = 0.03;
pos7 = pos6;
pos7(1) = pos6(1) + pos6(3) + pdix - expand;
pos7(3) = 1 - pdoy - shift(1) - (pos6(1) + pos6(3) + pdiy) + expand;
pos7(4) = 3*pos6(4) + 2*pdiy;
ax(7) = axes(figure(42),'position',pos7);

% boxes on
for i = 1:length(ax)
    box(ax(i),'on')
end

% get rid of interior time axes
for i = [1 2 4 5]
    %set(ax(i),'XTickLabel',{})
end

linkaxes(ax(1:6),'x')
xlim(ax(1),extrema(edges.time))

