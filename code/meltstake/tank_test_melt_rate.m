% Script to process ADV deployments in the glider tank 12 Jan 2024

clear

angle = [0 10 20];
pck(length(angle)) = struct;

data_path = '../../../../instrumentation/adv/data/pck_tank_test';

% load data
for i = 1:3
    load(fullfile(data_path,sprintf('pck%02d',angle(i)),'adv.mat'),'adv');
    pck(i).pck_time = adv.pck_time;
    pck(i).pck_dist = adv.pck_dist;
    pck(i).pck_amp = adv.pck_amp;
    pck(i).angle = angle(i);
    pck(i).pck_amp_plt = adv.pck_amp;
end

%% trim time frame to use same axis across all tests
tlim = [45 880];
for i = 1:3
    t_sec = seconds(pck(i).pck_time-pck(i).pck_time(1));
    idxt = t_sec >= tlim(1) & t_sec <= tlim(2);
    pck(i).pck_time = pck(i).pck_time(idxt);
    pck(i).pck_dist = pck(i).pck_dist(idxt,:);
    pck(i).pck_amp = pck(i).pck_amp(idxt,:,:);
    pck(i).pck_amp_plt = pck(i).pck_amp_plt(idxt,:,:);
end

%% distance from transducer to wall
m = 2.5/10./cosd(angle); % known melt rate (mm/s)
nt = 90;
dt = 10;
time_wall = dt*(0:nt-1)';
d0 = [80 55 55]; % need another offset due to transducer rotating toward wall

for i = 1:3
    pck(i).time_wall = time_wall;
    pck(i).dist_wall = (d0(i) + m(i)*time_wall)/cosd(angle(i));
end

%% clean up pck amplitude for better edge-finding
% ignore double reflections
r0 = 100;
r_max = r0 + 1.5*m(end)*seconds(pck(1).pck_time-pck(1).pck_time(1));

% clean up
for i = 1:3
    for j = 1:length(pck(i).pck_time)
        idxr = pck(i).pck_dist(j,:) > r_max(j);
        pck(i).pck_amp(j,idxr,:) = 0;
    end
end

%% calculate distance based on pck profiles
x0 = [50 50 50];
L = 100;
L0 = 50;
hann_width = 5;
for i = 1:3
    out = findADVIceEdge(pck(i),0,x0(i),L,L0,hann_width,30);
    pck(i).calc_time = out.time;
    pck(i).calc_pos = out.pos;
    pck(i).calc_desc = out.desc;
    pck(i).calc_hann_width = out.hann_width;
end

%% sigma filter edge position
nsigma = [2 2 2;2 2 2;2 1.9 2];
npass = [0 0 0;0 1 0;1 3 0];

for i = 1:3
    pck(i).calc_pos_filt = pck(i).calc_pos;
    for j = 1:3
        pos = pck(i).calc_pos(:,j);
        [~,idxs] = sigmaFilter(detrend(pos),nsigma(i,j),1,npass(i,j));
        pos(idxs) = nan;
        pck(i).calc_pos_filt(:,j) = pos;
    end
end

% manual trim beam2 (20deg)
pos = pck(3).calc_pos(:,2);
pos(pos<79) = nan;
pck(3).calc_pos_filt(:,2) = pos;


%% melt rate
m_obs = nan(3);
m_ci = nan(3);
r_fit = nan(3);

for i = 1:3
    tfit = seconds(pck(i).calc_time-pck(i).calc_time(1));
    for j = 1:3
        [b,bint,r,rint,stats] = regress(pck(i).calc_pos_filt(:,j),[ones(size(tfit)) tfit],0.05);
        m_obs(i,j) = round(b(2),3);
        m_ci(i,j) = round(diff(bint(2,:))/2,3);
        r_fit(i,j) = round(b(1),1);
    end
end

%% calculate error
m_obs_mean = mean(m_obs,2)';
error = (m_obs_mean-m)./m;

%% plot
fs  = 11;
lw = 1.5;

figsize = [16 16];
fig = figure(6); clf
setFigureSize(fig,figsize);
set(fig,'inverthardcopy','off','color',[1 1 1])

% limits
CLIM = [55 200];
RLIM = [3 35];

% axes
clear ax
pad = [.02 .04 .14 .08];
shift = [0.02 -0.01];

bm_lbl = {'b1 (down)','b2 (near)','b3 (far)'};

clear ax
for i = 1:3
    for j = 1:3
        ax(3*(i-1)+j) = axes(figure(6),'position',axgridpos(3,3,3*(i-1)+j,pad,shift));
        box on
        hold on
        
        % pcolor
        t0 = pck(i).pck_time(1);
        t = seconds(pck(i).pck_time-t0);
        tmax = max(t);
        pcolor(t,pck(i).pck_dist(2,:)/10,pck(i).pck_amp_plt(:,:,j)')
        shading flat
        
        % calculated
        t = seconds(pck(i).calc_time-t0);
        plot(t,pck(i).calc_pos(:,j)/10,'r-','linewidth',lw)
        plot(t,pck(i).calc_pos_filt(:,j)/10,'k-','linewidth',lw+0.5)

        % fit
        plot(t,(r_fit(i,j)+m_obs(i,j)*t)/10,'w--','linewidth',lw)
        text(500,6,sprintf('%.3f%c%.3f mm/s',round(m_obs(i,j),3),char(177),round(m_ci(i,j),3)),...
            'fontsize',fs-2,'color','w','horizontalalignment','center')
        
        % actual
        t = pck(i).time_wall;
        plot(t,pck(i).dist_wall/10+5,'-.','color','w','linewidth',lw-0.2)

        % panel labels
        text(ax(3*(i-1)+j),.01,.98,sprintf('(%c%d)',char(96+i),j),'units','normalized',...
            'fontsize',fs-1,'verticalalignment','top','horizontalalignment','left',...
            'color','w')
    end
end

%%% general %%%
linkaxes(ax)
clim(CLIM)
ylim(RLIM)

% clean up axes labels
for i = 1:6
    set(ax(i),'xticklabel',{})
end
for i = setxor(1:9,1:3:9)
    set(ax(i),'yticklabel',{})
end
for i = 1:3
    ylabel(ax(3*(i-1)+1),...
        {sprintf('angle: %d^\\circ',angle(i)),sprintf('m_{true}: %.3f mm/s',m(i)),'---------------------------','ADV distance [cm]'},'fontsize',fs-1)
end
for i = 7:9
    xlabel(ax(i),'time [s]','fontsize',fs)
end

% column titles
for i = 1:3
    title(ax(i),sprintf('ADV receiver %d',i),'fontsize',fs)
end

% set global x limit
xlim(ax(1),[0 diff(tlim)])

% colorbar
cbar = colorbar(ax(3),'position',cbarpos(ax(3),.01,.02));
cbar.Label.String = 'amplitude [counts]';
cbar.Label.FontSize = fs;

% legend
h = findobj(ax(1),'type','line');
h = h([3 4 2 1]);
lgd_lbls = {'edge (kept)','edge (rejected)','fit','true melt rate'};
lgd = legend(ax(3),h,lgd_lbls,'fontsize',fs-2,'color',0.7*[1 1 1],...
    'orientation','horizontal');
lgd.Position = [0.165 0.95 0.7050 0.0331];
        
        