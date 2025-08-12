% Script to create figure 3 for the melt rate GRL paper. This figure is a
% multi-panel scatter plot comparing melt rate observations to melt rates
% predicted by various parameterizations. Output for parameterizations is
% created using the script "ms_melt_models.m"
%
% KJW
% 7 Apr 2024
set(0,'defaulttextinterpreter','latex');

clear

% load data
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

% sections
sections = [8 3 29 9];
lbls = {'A','B','C','D'};
mkrs = {'>','^','square','hexagram'};
lgd_lbls = {'A (crossflow)','B (plume)','C (wave)','D (eddy)'};

% observations
m_obs = msTable2Vector(ms_tbl.m);
m_ci = ms_tbl.m_ci;
m_obs = m_obs*.24; % cm/hr --> m/day
m_ci = m_ci*.24;
[u,u_ci] = msTable2Vector(ms_tbl.u0);
up = ms_tbl.uprimesqrd;
[T,T_ci] = msTable2Vector(ms_tbl.T);

m_obs(end) = nan;
m_ci(end) = nan;

% models
load ../../data/melt_models.mat

% melt error
r_3eqn = m_obs./m_mean(:,1);

r = repmat(m_obs,[1 size(m_mean,2)])./m_mean;
r_ci = repmat(m_ci,[1 size(m_mean,2)])./m_mean;
%% plot
fs  = 11;
lw = 1;
clear ax
clear h

% padding
pxi = .05;
pxo = .15;
pyi = .05;
pyo = .15;

% create and size figure
figsize = [9 10];
figure(3); clf
set(figure(3),'units','centimeters'); 
fpos = get(figure(3),'position');
fpos(3:4) = figsize;
set(figure(3),'position',fpos);
set(figure(3),'paperunits','centimeters')
set(figure(3),'papersize',figsize)

panel_lbls = {'3-equation','convective','modified 3-eqn','bulk iceberg'};

for i = 1:1
    % axis
    ax(i) = axes(figure(3),'position',axgridpos(1,1,i,pxi,pyi,pxo,pyo));
    hold on
    box on

    % plot
    plot([0 15],1*[0 15],'k-')
    plot([0 15],2*[0 15],'k--')
    plot([0 15],4*[0 15],'k--')
    errorbar(m_mean(:,i),m_obs,-m_ci,m_ci,...
        'color',.5*[1 1 1],'linestyle','none','linewidth',lw)
    plot(m_mean(:,i),m_obs,'ko','markersize',4,'markerfacecolor',colors(1))
    xlim([0 2])
    ylim([0 2])
%     for j = 1:length(sections)
%         plot(m_mean(sections(j),i),m_obs(sections(j)),'ko','markersize',4,'markerfacecolor','r')
%         text(ax(i),m_mean(sections(j),i),m_obs(sections(j)),lbls{j},'fontsize',fs)
%     end


    if any(i==[1 3])
        ylabel('observed melt rate (m/day)','fontsize',fs)
    end
    if any(i==[1 3 4])
        xlabel('modeled melt rate (m/day)','fontsize',fs)
    end

%     text(ax(i),.01,.98,['(' char(96+i) ') ' panel_lbls{i}],'units','normalized','fontsize',fs,...
%         'verticalalignment','top','horizontalalignment','left')
end

% add sections
for i = 1:1
    for j = 1:length(sections)
        h(j) = plot(ax(i),m_mean(sections(j),i),m_obs(sections(j)),'ko','markersize',7,...
            'markerfacecolor',colors(j+1),'marker',mkrs{j});
    end
end
legend(ax(1),h,lgd_lbls,'location','southeast','fontsize',fs-2,'color','w')

true_aspect

% xlim(ax(3),[0 15])
% ylim(ax(3),[0 15])

