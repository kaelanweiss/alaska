% Script to create methods/data examples figure for glacier GRL.
%
% KJW
% 24 Sep 2025
clear

addpath('..')
load glacier_clrs.mat
load F:/adv/mean_profile.mat
proc_dir = 'F:/meltstake/data/proc';
raw_dir = 'F:/meltstake/data/raw';

ms_tbl = loadMSInfo(26:28,'segments');
dep_nums = unique(ms_tbl.Number);

% load segments
segments = [4 6 2];
ms(length(segments)) = struct();
for i = 1:length(segments)
    dep_name = ms_tbl.Folder{find(ms_tbl.Number==(i+25),1)};
    row_idx = ms_tbl.Number==(i+25) & ms_tbl.Window==segments(i);
    % things
    ms(i).rmax = ms_tbl.rmax(row_idx);
    ms(i).depth = ms_tbl.depth(row_idx);
    ms(i).beam = str2double(strsplit(ms_tbl.Beam{row_idx},','));
    ms(i).mb = [msTable2Vector(ms_tbl.m1(row_idx)) msTable2Vector(ms_tbl.m2(row_idx)) msTable2Vector(ms_tbl.m3(row_idx))]*0.24;
    ms(i).m = msTable2Vector(ms_tbl.m(row_idx))*0.24;
    [ci1,ci2] = msTable2Vector(ms_tbl.m_ci(row_idx));
    ms(i).m_ci = [ci1 ci2]*0.24;
    % adcp
    dat = load(fullfile(proc_dir,dep_name,sprintf('adcp%d.mat',segments(i))));
    ms(i).adcp = dat.adcp;
    % T
    dat = load(fullfile(proc_dir,dep_name,sprintf('T%d.mat',segments(i))));
    ms(i).T = dat.T;
    % adv
    dat = load(fullfile(raw_dir,dep_name,'adv','adv.mat'),'adv');
    idxt = dat.adv.pck_time>=ms_tbl.Start(row_idx) & dat.adv.pck_time<=ms_tbl.End(row_idx);
    ms(i).adv = struct('pck_time',dat.adv.pck_time(idxt),'pck_dist',dat.adv.pck_dist(idxt,:),'pck_amp',dat.adv.pck_amp(idxt,:,:));
    
    ms(i).adv.pck_amp_norm = ms(i).adv.pck_amp./repmat(mean_prof',[length(ms(i).adv.pck_time) 1 3]);
    ms(i).adv.pck_amp_norm = 100*ms(i).adv.pck_amp_norm/max(ms(i).adv.pck_amp_norm,[],'all');
end

%% plot
lw = 1;
fs = 11;

fig_size = [14 12];
fig = figure(1); clf
setFigureSize(fig,fig_size);
clear ax

pad = [.02 .04 .12 .07];
shift = [0 0.02];

% limits
VLIM = 0.4;
TLIM = [3.5 7.5];
RMAX = 0.55;
PLIM = [0.1 0.25];
BLIM = {[110 210],[100 180],[125 210]};

% panel labels
pnl_lbls = {'u (horizontal, along-ice)','w (vertical)','ocean temperature','ADV backscatter'};

% plot melt rate
d0 = [0.17 0.2 0.16];

adv_bn = [1 3 2];

% create axes
for i = 1:12
    ax(i) = axes(fig,'position',axgridpos(4,3,i,pad,shift,'flip'),'fontsize',fs-2);
    box(ax(i),'on')
    hold(ax(i),'on')
end

% pplort
for i = 1:3
    ax0 = 4*(i-1)+1;
    % u
    y = ms(i).rmax/cosd(25) - ms(i).adcp.burst.range;
    pcolor(ax(ax0),ms(i).adcp.burst.time,y,ms(i).adcp.burst.vel_ice(:,:,1)')
    shading(ax(ax0),'flat')
    cmocean('bal',ax(ax0))

    % w
    pcolor(ax(ax0+1),ms(i).adcp.burst.time,y,ms(i).adcp.burst.vel_ice(:,:,2)')
    shading(ax(ax0+1),'flat')
    cmocean('bal',ax(ax0+1))

    % T
    for j = 1:3
        h(j) = plot(ax(ax0+2),ms(i).T.time,ms(i).T.T(:,j),'color',T_clrs(4-j,:),'linewidth',lw);
    end

    % ice
    pcolor(ax(ax0+3),ms(i).adv.pck_time,ms(i).adv.pck_dist(2,:)/1000,ms(i).adv.pck_amp(:,:,adv_bn(i))')
    shading(ax(ax0+3),'flat')
%     cmocean('ice',ax(ax0+3))
    dt = diff(ms(i).adcp.burst.time([1 end]));
    plot(ax(ax0+3),[0 dt]+ms(i).adcp.burst.time(1),d0(i)+[0 days(dt)]*ms(i).m,'k--','linewidth',lw)
    text(ax(ax0+3),ms(i).adcp.burst.time(1)+dt/2,PLIM(1),sprintf('%.1f [%.1f,%.1f] m/dy',round(ms(i).m,1),round(ms(i).m-ms(i).m_ci(1),1),round(ms(i).m+ms(i).m_ci(2),1)),...
        'verticalalignment','bottom','horizontalalignment','center','fontsize',fs-2)
    linkaxes(ax(ax0:(ax0+3)),'x')

    % panel label
    for j = 1:4
        text(ax(ax0+j-1),0.015,1,sprintf('(%c)',char(96+4*(i-1)+j)),'units','normalized','fontsize',fs-1,...
            'verticalalignment','top','horizontalalignment','left')
    end

    % row label
    if i==2
        for j = 1:4
            text(ax(ax0+j-1),0.5,1.01,pnl_lbls{j},'units','normalized','fontsize',fs,...
                'verticalalignment','bottom','horizontalalignment','center')
        end
    end
end

% set axes limits
% time limits
for i = 1:3
    xlim(ax(4*(i-1)+1),extrema(ms(i).adcp.burst.time))
end
% adcp limits
for i = [1 2 5 6 9 10]
     ylim(ax(i),[0 RMAX])
     clim(ax(i),VLIM*[-1 1])
     set(ax(i),'ytick',0:0.25:1)
end
% T limits
for i = [3 7 11]
    ylim(ax(i),TLIM)
    set(ax(i),'ytick',3:10)
end
% adv
for i = [4 8 12]
    ylim(ax(i),PLIM)
    clim(ax(i),BLIM{i/4})
end

% clear unwanted labels
for i = setxor(1:12,[4 8 12])
    set(ax(i),'xticklabels',{})
end
for i = 5:12
    set(ax(i),'yticklabels',{})
end

% labels
for i = [1 2]
    ylabel(ax(i),'y [m]','fontsize',fs)
end
ylabel(ax(3),'T_w [^\circC]','fontsize',fs)
ylabel(ax(4),'ADV distance [m]','fontsize',fs)

% T legend
lgd = legend(ax(11),flip(h),{'0.40 m','0.25 m','0.10 m'},'location','southeast','fontsize',fs-2);
lgd.Position(1:2) = [0.8 0.3];

% velocity colorbar
cbar = colorbar(ax(9),'position',cbarpos(ax(9:10),.01,.025));
text(ax(9),1.17,1.16,'m/s','fontsize',fs-2,'units','normalized')

% TO DO
% ice patch in adcp plots


