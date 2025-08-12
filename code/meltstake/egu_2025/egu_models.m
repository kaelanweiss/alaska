% Script to plot melt rate model results for EGU 2025
%
% KJW

clear

addpath('..')

% load
load clrs_egu.mat
mod_glac = load('../../../data/melt_models_glacier.mat');
mod_berg = load('../../../data/melt_models.mat');
mod_tau = load('melt_3eqn_tau.mat');
ms_tbl = loadMSInfo('manualwindows');

% separate glacier and berg data
iglc = ms_tbl.Number >= 26;
ibrg = ~iglc;

% separate individual glacier deployments
idx_dep = false([length(iglc),3]);
for i = 1:3
    idx_dep(:,i) = ms_tbl.Number == i+25;
end

% observed melt
m_obs = msTable2Vector(ms_tbl.m);
[m_ci_l,m_ci_h] = msTable2Vector(ms_tbl.m_ci);
m_obs = m_obs*0.24; % cm/hr --> m/day
m_ci_l = m_ci_l*0.24;
m_ci_h = m_ci_h*0.24;

% calculate statistics for error vs sample tau
r = repmat(m_obs(iglc),[1 size(mod_tau.m_mean,[2 3])])./mod_tau.m_mean;
r_med = squeeze(median(r,'omitnan'));
r_mean = squeeze(mean(r,'omitnan'));
r_std = squeeze(std(r,'omitnan'));
r_fit = nan(size(r_med));
r_fit_ci = nan(size(r_med));
tau_lbls = cell(size(mod_tau.tau));
for i = 1:size(r,2)
    for j = 1:2
        [r_fit(i,j),dr] = regress(m_obs(iglc),mod_tau.m_mean(:,i,j));
        r_fit_ci(i,j) = diff(dr)/2;
    end
    if i<=length(mod_tau.tau)
        tau_lbls{i} = sprintf('2^{%d}',log2(mod_tau.tau(i)));
    end
end

%% plot
fs = 11;
lw = 1;
ms = 7;
fig1 = setFigureSize(figure(1),[9 8.5]);
fig2 = setFigureSize(figure(2),[9 8.5]);
fig3 = setFigureSize(figure(3),[9 8.5]);
clf(fig1); clf(fig2); clf(fig3)
figs = [fig1 fig2 fig3];
smbs = {'o','square','^'};
clr_brg = 0.6*[1 1 1];

% m_obs vs m_mod scatter plot
mod_nums = [1 5 6];
mod_ttls = {'full resolution','coarse resolution','modified (Schulz et al., 2022)'};
clear ax
for i = 1:3
    clear h
    ax(i) = axes(figs(i));
    hold(ax(i),'on')
    box(ax(i),'on')
    plot(ax(i),[0 20],[0 20],'k-')
    % error bars
    errorbar(ax(i),[mod_berg.m_mean(:,mod_nums(i)); mod_glac.m_mean(:,mod_nums(i))],m_obs,-m_ci_l,m_ci_h,'color',.5*[1 1 1],'linestyle','none','linewidth',lw)
    % 3eqn (full)
    h(4) = plot(ax(i),mod_berg.m_mean(:,mod_nums(i)),m_obs(ibrg),'ko','markerfacecolor',clr_brg,'markersize',ms-3);
    for j = 1:3
        h(j) = plot(ax(i),mod_glac.m_mean(idx_dep(iglc,j),mod_nums(i)),m_obs(idx_dep(:,j)),['k' smbs{j}],'markerfacecolor',clrs(j,:),'markersize',ms);
    end
    % ensemble metrics
%     if i<3
%         text(ax(i),.7,.55,{sprintf('r_{med}=%.1f',round(r_med(1,i),1)),sprintf('r_{fit}=%.1f%c%.1f',round(r_reg(1,i),1),char(177),round(diff(r_reg_int)/2,1))},'units','normalized','fontsize',fs)
%     else
%         [iint_i] = regress(m_obs(iglc),mod_glac.m_mean(:,mod_nums(i)));
%         r_med_i = median(m_obs(iglc)./mod_glac.m_mean(:,mod_nums(i)),'omitnan');
%         text(ax(i),.7,.55,{sprintf('r_{med}=%.2f',round(r_med,2)),sprintf('r_{fit}=%.2f%c%.2f',round(r_reg,2),char(177),round(diff(r_reg_int)/2,2))},'units','normalized','fontsize',fs)
%     end
    % labels
    xlabel(ax(i),'m_{3eqn} [m/day]','fontsize',fs)
    ylabel(ax(i),'m_{obs} [m/day]','fontsize',fs)
    title(ax(i),mod_ttls{i},'fontsize',fs+2)
    xlim(ax(i),[0 4.3])
    ylim(ax(i),[0 4.3])
    
    % histogram
    r_hist = m_obs(iglc)./mod_glac.m_mean(:,mod_nums(i));
    ax_in = axes(figs(i),'position',[0.645 0.638-.05 0.24 0.26],'units','normalized','fontsize',fs-2);
    hold(ax_in,'on')
    box(ax_in,'on')
%     histogram(ax_in,log10([r_hist; m_obs(ibrg)./mod_berg.m_mean(:,mod_nums(i))]),-1.5:0.15:2.3,'facecolor',0*colors(7),'facealpha',1)
    histogram(ax_in,log10(r_hist),-1.5:0.15:2.3,'facecolor',0.*colors(3),'facealpha',1)
    plot(ax_in,log10(median(r_hist,'omitnan'))*[1 1],[0 15],'k-','linewidth',1)
    xlim(ax_in,[-1 2])
    ylim(ax_in,[0 12])
    set(ax_in,'xtick',-2:2)
    set(ax_in,'xticklabelrotation',0)
    grid(ax_in,'on')

    xlabel(ax_in,'$\textsf{log}_\textsf{10}\Big(\frac{\textsf{m}_\textsf{obs}}{\textsf{m}_\textsf{3eqn}}\Big)$','fontsize',fs-2,'interpreter','latex')
%     ylabel(ax_in,'occur.','fontsize',fs)

    if i == 1
        legend(ax(1),h,{'dep 1','dep 2','dep 3','icebergs'},'fontsize',fs-1,'location','southeast')
    end
end
xlim(ax(3),[0 12])
ylim(ax(3),[0 12])

%% r vs tau
ms = 7;
fig4 = setFigureSize(figure(4),[11 7.5]);
clf(fig4)
clear h2
ax2 = axes(fig4);
hold(ax2,'on')
box(ax2,'on')
h2(1) = plot(mod_tau.tau,r_mean(1:end-1,1),'ko','markerfacecolor','k','markersize',ms);
h2(2) = plot(mod_tau.tau,r_mean(1:end-1,2),'ksquare','markerfacecolor',0.5*[1 1 1],'markersize',ms+1);
set(ax2,'xscale','log')
set(ax2,'xtick',10.^(-1:1:4))
xlim(ax2,extrema(mod_tau.tau).*[0.5 2])
ylim(ax2,[1 9])
% set(ax2,'xtick',log2(mod_tau.tau))
% set(ax2,'xticklabels',tau_lbls)
xlabel(ax2,'down-sampled forcing resolution [s]','fontsize',fs)
ylabel(ax2,'mean( m_{obs}/m_{3eqn} )','fontsize',fs)
legend(ax2,h2,{'free stream','mean components'},'fontsize',fs,'location','northwest')



