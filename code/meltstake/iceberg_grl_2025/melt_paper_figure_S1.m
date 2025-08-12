% Script to create a figure showing the process for determining melt rate
% from ADV profiles.
%
% Rows correspond to:
%   1. raw backscatter amplitude
%   2. backscatter + edgefinding
%   3. 2 + filter
%   4. 3 + linear fit
%
% KJW
% 29 Dec 2024

clear

addpath('..')

% set deployment number
dep_num = 8;

% meltstake table
ms_tbl = loadMSInfo('manualwindows');
ms_filt = loadMSInfo('melt_filter');

% get deployment name, number of segments, segment times
dep_name = ms_tbl.Folder{find(dep_num==ms_tbl.Number,1)};
n_seg = sum(ms_tbl.Number==dep_num);
t_start = ms_tbl.Start(ms_tbl.Number==dep_num);
t_end = ms_tbl.End(ms_tbl.Number==dep_num);

% filter params
ms_filt = ms_filt(ms_filt.Number==dep_num,:);
nsigma = [ms_filt.n_sigma1 ms_filt.n_sigma2 ms_filt.n_sigma3];
npass = [ms_filt.n_pass1 ms_filt.n_pass2 ms_filt.n_pass3];

% melt rates
m1 = ms_tbl.m1(ms_tbl.Number==dep_num);
m2 = ms_tbl.m2(ms_tbl.Number==dep_num);
m3 = ms_tbl.m3(ms_tbl.Number==dep_num);
m = [msTable2Vector(m1) msTable2Vector(m2) msTable2Vector(m3)];

% load data
raw_path = fullfile('F:meltstake\data\raw',dep_name);
proc_path = fullfile('F:meltstake\data\proc',dep_name);
load(fullfile(raw_path,'adv','adv.mat'),'adv')
load(fullfile(proc_path,'adv','pck_edges.mat'))

%% plot
fs  = 11;
lw = 1;

figsize = [16 18];
fig = figure(5); clf
setFigureSize(fig,figsize);
set(fig,'inverthardcopy','off','color',[1 1 1])

% axes
clear ax
pad = [.01 .04 .09 .04];
shift = [-0.03 0.01];
nrows = 5;
npans = nrows*3;

% limits
CLIM = [55 200];
RLIM = [8 30];

% panel labels
panel_lbls = {'echosounder','edge-finding','drill correction','filter','linear fit'};

% create axes and plot pcolor
for i = 1:3
    for j = 1:nrows
        % axis
        ax(i+3*(j-1)) = axes(fig,'position',axgridpos(nrows,3,i+3*(j-1),pad,shift));
        % pcolor
        pcolor(ax(i+3*(j-1)),adv.pck_time,adv.pck_dist(1,:)/10,adv.pck_amp(:,:,i)')
        shading flat
        hold(ax(i+3*(j-1)),'on');
        clim(CLIM)
        % segments
        for k = 1:n_seg
            plot(ax(i+3*(j-1)),[t_start(k) t_end(k)]+minutes(1)*[1 -1],RLIM(1)*[1 1]+1,'|-',...
            'linewidth',lw-0.25,'markersize',5,'color',0.8*[1 1 1])
            if j == 1
                text(mean([t_start(k) t_end(k)]),RLIM(1)+1,sprintf('seg%d',k+find(ms_tbl.Number==dep_num,1)-1),...
                    'verticalalignment','bottom','horizontalalignment','center',...
                    'fontsize',fs-2,'color',1.*[1 1 1],'fontweight','bold')
            end
        end
        % panel labels
        text(ax(i+3*(j-1)),.01,.98,sprintf('(%c%d) %s',char(96+j),i,panel_lbls{j}),'units','normalized',...
            'fontsize',fs-1,'verticalalignment','top','horizontalalignment','left')
    end
end

linkaxes(ax)
ylim(ax(1),RLIM)

%%% row 2 %%%
for i = 1:3
    plot(ax(i+3),edges.time,edges.pos(:,i)/10,'k-','linewidth',lw)
end

%%% row 3 %%%
% fix jump
jump = [87.4 89.1 87.4];
edges.pos_jump = edges.pos;
idxt = edges.time <= tdrill;
for i = 1:3
    edges.pos_jump(idxt,i) = edges.pos_jump(idxt,i) - jump(i);
end
% plot
for i = 1:3
    plot(ax(i+6),edges.time,edges.pos_jump(:,i)/10,'k-','linewidth',lw)
end

%%% row 4 %%%
edges.pos_filt = nan*edges.pos_jump;
% filter segments
for j = 1:n_seg
    idxt = edges.time >= t_start(j) & edges.time < t_end(j);
    for i = 1:3
        pos = edges.pos_jump(idxt,i);
        if npass(j,i) ~= -1
            [~,idxs] = sigmaFilter(detrend(pos),nsigma(j,i),1,npass(j,i));
            pos(idxs) = nan;
            edges.pos_filt(idxt,i) = pos;
        end
    end
end
% plot
for i = 1:3
    if i == 1
        h4(1) = plot(ax(i+9),edges.time,edges.pos_jump(:,i)/10,'r-','linewidth',lw);
        h4(2) = plot(ax(i+9),edges.time,edges.pos_filt(:,i)/10,'k-','linewidth',lw);
    else
        plot(ax(i+9),edges.time,edges.pos_jump(:,i)/10,'r-','linewidth',lw)
        plot(ax(i+9),edges.time,edges.pos_filt(:,i)/10,'k-','linewidth',lw)
    end
end
% legend
lgd4 = legend(ax(12),h4,{'reject','keep'},'fontsize',fs-2,'units','normalized',...
    'position',[0.8 ax(12).Position(2)+0.13 0.12 0.02],'color',0.8*[1 1 1]);

%%% row 5 %%%
for i = 1:3
    plot(ax(i+12),edges.time,edges.pos_filt(:,i)/10,'k-','linewidth',lw+0.5)
    for j = 1:n_seg
        idxt = edges.time >= t_start(j) & edges.time <= t_end(j);
        dt = hours(t_end(j)-t_start(j));
        epst = .04*dt;
        if ~isnan(m(j,i))
            y_bar = mean(edges.pos_filt(idxt,i)/10,'omitnan');
            dy = m(j,i)*(dt-2*epst);
            plot(ax(i+12),[t_start(j) t_end(j)]+hours(epst)*[1 -1],dy*[-1 1]/2+y_bar-0,'w-','linewidth',lw)
            text(ax(i+12),mean([t_start(j) t_end(j)]),y_bar-3,sprintf('%.2f m/d',m(j,i)*0.24),...
                'color','w','horizontalalignment','center','fontsize',fs-2)
        end
    end
end
% legend
h5 = cat(1,findobj(ax(15),'type','line','color','k'),findobj(ax(15),'type','line','color','w'));
lgd5 = legend(ax(15),h5(1:2),{'keep','fit'},'fontsize',fs-2,'units','normalized',...
    'position',[0.8 ax(15).Position(2)+0.13 0.12 0.02],'color',0.8*[1 1 1]);

%%% general %%%
% clean up axes labels
for i = 1:3*(nrows-1)
    set(ax(i),'xticklabel',{})
end
for i = setxor(1:npans,1:3:npans)
    set(ax(i),'yticklabel',{})
end
for i = 1:3:npans
    ylabel(ax(i),'ADV distance [cm]','fontsize',fs-1)
end
for i = npans-2:npans
    set(ax(i),'xticklabelrotation',0)
    xtickformat(ax(i),'HH:mm')
end

% column titles
for i = 1:3
    title(ax(i),sprintf('ADV receiver %d',i),'fontsize',fs)
end

% set global x limit
xlim(ax(1),[t_start(1) t_end(end)])

% colorbar
cbar = colorbar(ax(3),'position',cbarpos(ax(3),.01,.02));
cbar.Label.String = 'amplitude [counts]';
cbar.Label.FontSize = fs;