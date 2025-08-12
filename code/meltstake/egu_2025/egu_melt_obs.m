% Script to create melt observations figure for EGU 2025
%
% KJW
% 7 Apr 2025
%
% USE SYMBOLS TO DIFFERENTIATE DEPLOYMENTS
set(0,'defaulttextinterpreter','latex');

clear

addpath('..')
load clrs_egu

% load data
ms_tbl = loadMSInfo(26:28,'manualwindows');
ms_tbl = ms_tbl(1:end-1,:);

dep_nums = unique(ms_tbl.Number);

% separate individual glacier deployments
idx_dep = false([size(ms_tbl,1),length(dep_nums)]);
for i = 1:3
    idx_dep(:,i) = ms_tbl.Number == i+25;
end

% cm/hr to m/dy
cmh2md = 0.24;

% extract fields
[m,~] = msTable2Vector(ms_tbl.m);
[T,T_ci] = msTable2Vector(ms_tbl.T);
[S,S_ci] = msTable2Vector(ms_tbl.S);
[u,u_ci] = msTable2Vector(ms_tbl.u_max);
[m_ci_l,m_ci_h] = msTable2Vector(ms_tbl.m_ci);

% convert melt rate to m/day
m = m*cmh2md;
m_ci_l = m_ci_l*cmh2md;
m_ci_h = m_ci_h*cmh2md;

% standard error
tau_decorr = ms_tbl.tau_decorr;
Neff = 60*ms_tbl.Duration./tau_decorr;
NS = ms_tbl.nS;

T_std = T_ci;
u_std = u_ci;
S_std = S_ci;

T_ci = T_ci./sqrt(Neff);
S_ci = S_ci./sqrt(NS);
u_ci = u_ci./sqrt(Neff);

%% UT contours
n_cont = 200;
u_cont = linspace(0,.3,n_cont);
t_cont = linspace(0,11,n_cont);
[U_cont,T_cont] = meshgrid(u_cont,t_cont);
UT_cont = U_cont.*T_cont;
m_cont = solve3Eqn(U_cont,T_cont,25,0)*86400;
cont_lvls = .05:.2:5;

%% plot
lw = 1;
fs = 11;
marker_scale = 75;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2.1 parameter space (u vs T with colors corresponding to S)
figure(21); clf
% create and size figure
figsize = [10 10];
setFigureSize(figure(21),figsize);

% plot
clear h
ax1 = axes(figure(21),'fontsize',fs-2);
hold on
contour(u_cont,t_cont,m_cont,cont_lvls,'k:')
errorbar(u,T,-T_ci,T_ci,-u_ci,u_ci,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
for i = 1:length(dep_nums)
    h(i) = scatter(u(idx_dep(:,i)),T(idx_dep(:,i)),marker_scale,m(idx_dep(:,i)),smbs{i},'filled','markeredgecolor','k');
end

cbar = colorbar;
cbar.Location = 'northoutside';
box on
legend(ax1,h,{'dep 1','dep 2','dep 3'},'location','southwest','fontsize',fs-1)

xlim([.05 .3])
ylim([5 7])
clim([1.5 4])

% labels
xlabel('$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','fontsize',fs)
ylabel('$\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}$ \textsf{[K]}','fontsize',fs)
cbar.Label.String = 'm_{obs} [m/day]';
cbar.Label.FontSize = fs;
cmocean('-ice')

% contour text label
text(0.26,5.6,'$\overline{\textsf{U}_\textsf{max}}\cdot\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}\textsf{=}$\textsf{const.}','rotation',-72,...
    'verticalalignment','middle','horizontalalignment','center')

% tighten up colorbar spacing
axpos = ax1.Position;
ax1.Position = axpos + [0 0 0 .05];
cbar.Position(4) = .04;
cbar.Position(2) = sum(ax1.Position([2 4])) + .012;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 2.2 (m vs u,T,S)
figure(22); clf
% create and size figure
figsize = [9 11];
setFigureSize(figure(22),figsize);

x_data = {u,T,S};
x_ci = {u_ci,T_ci,S_ci};
x_lbls = {'$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','$\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}$ \textsf{[K]}','$\overline{\textsf{S}_\textsf{w}}$ \textsf{[psu]}'};
c_data = [2,1,2];
c_cmap = {'-solar','tempo','-solar'};
c_lbls = {'$\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}$ \textsf{[K]}','$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','$\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}$ \textsf{[K]}'};
c_lims = {[0.05 .3],[5 7]};
clear ax2
for i = 1:3
    ax2(i) = axes(figure(22),'position',axgridpos(3,1,i,[.1,.12,.15,.1],[-.04 0]),'fontsize',fs-2);
    hold on
    errorbar(x_data{i},m,-m_ci_l,m_ci_h,-x_ci{i},x_ci{i},'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
    for j = 1:length(dep_nums)
        scatter(x_data{i}(idx_dep(:,j)),m(idx_dep(:,j)),18,x_data{c_data(i)}(idx_dep(:,j)),'filled',smbs{j},'markeredgecolor','k')
    end
    xlabel(x_lbls{i},'fontsize',fs)
%     ylabel('$\textsf{m}_\textsf{obs}$ \textsf{[m/day]}','fontsize',fs)
    ylabel('m_{obs} [m/day]','fontsize',fs,'interpreter','tex')
    ylim([0 4.5])
    box on
    clim(c_lims{c_data(i)})
    if i==3
        xlim([22 29])
        set(gca,'clipping','off')
    end

    cbar = colorbar;
    cbar.Position = cbarpos(gca,.005,.02);
    cbar.Label.String = c_lbls{i};
    cbar.Label.FontSize = fs;
    cbar.Label.Interpreter = 'latex';
    cmocean(c_cmap{i})
end

xlim(ax2(1),[0.05 0.28])
xlim(ax2(2),[5 7])
xlim(ax2(3),[25.5 28])

% %% plot without color data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure 2.1 parameter space (u vs T with colors corresponding to S)
% figure(23); clf
% % create and size figure
% figsize = [10 10];
% setFigureSize(figure(23),figsize);
% 
% % plot
% ax1 = axes(figure(23),'fontsize',fs-2);
% hold on
% contour(u_cont,t_cont,m_cont,cont_lvls,'k:')
% errorbar(u,T,-T_ci,T_ci,-u_ci,u_ci,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
% for i = 1:length(dep_nums)
%     idx = ms_tbl.Number == dep_nums(i);
%     plot(u(idx),T(idx),['k' smbs{i}],'markerfacecolor',clrs(i,:),'markersize',8)
% end
% box on
% 
% xlim([.05 .3])
% ylim([5 7])
% 
% % labels
% xlabel('$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','fontsize',fs)
% ylabel('$\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}$ \textsf{[K]}','fontsize',fs)
% 
% % contour text label
% text(0.26,5.6,'$\overline{\textsf{U}_\textsf{max}}\cdot\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}\textsf{=}$\textsf{const.}','rotation',-72,...
%     'verticalalignment','middle','horizontalalignment','center')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % figure 2.2 (m vs u,T,S)
% figure(24); clf
% % create and size figure
% figsize = [9 11];
% setFigureSize(figure(24),figsize);
% 
% x_data = {u,T,S};
% x_ci = {u_ci,T_ci,S_ci};
% x_lbls = {'$\overline{\textsf{U}_\textsf{max}}$ \textsf{[m/s]}','$\overline{\textsf{T}_\textsf{w}-\textsf{T}_\textsf{i}}$ \textsf{[K]}','$\overline{\textsf{S}_\textsf{w}}$ \textsf{[psu]}'};
% 
% clear ax2
% for i = 1:3
%     ax2(i) = axes(figure(24),'position',axgridpos(3,1,i,[.1,.12,.15,.1],[-.04 0]),'fontsize',fs-2);
%     hold on
%     errorbar(x_data{i},m,-m_ci_l,m_ci_h,-x_ci{i},x_ci{i},'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
%     for j = 1:length(dep_nums)
%         idx = ms_tbl.Number == dep_nums(j);
%         plot(x_data{i}(idx),m(idx),['k' smbs{j}],'markerfacecolor',clrs(j,:),'markersize',5)
%     end
%     xlabel(x_lbls{i},'fontsize',fs)
% %     ylabel('$\textsf{m}_\textsf{obs}$ \textsf{[m/day]}','fontsize',fs)
%     ylabel('m_{obs} [m/day]','fontsize',fs,'interpreter','tex')
%     ylim([0 4.5])
%     box on
% 
%     if i==3
%         xlim([22 29])
%         set(gca,'clipping','off')
%     end
% end
% 
% xlim(ax2(1),[0.05 0.28])
% xlim(ax2(2),[5 7])
% xlim(ax2(3),[25.5 28])