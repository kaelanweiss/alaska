% Script to plot CTD profiles of T, S, rho, and N for supplementaryu
% information.
%
% KJW
% 22 Oct 2026

clear

load F:/meltstake/data/proc/glacier_2024_ctd.mat
load glacier_clrs.mat

%% plot
fig = figure(1); clf
setFigureSize(fig,[12 6]);

fs = 10;
lw = 1.2;

pad = [.05 .05 .1 .1];
shift = [-0.02 0.07];
for i = 4:-1:1
    ax(i) = axes(fig,'position',axgridpos(1,4,i,pad,shift));
    hold on
end

% panel labels
x_lbls = {'T [\circC]','S [psu]','\rho [kg m^{-3}]','N/2\pi [cps]'};

% plot
% T
for i = 1:3
    plot(ax(1),ctd(i).T_mean,ctd(i).depth_mean,'-','color',dep_clrs(i,:),'linewidth',lw);
end
% S
for i = 1:3
    plot(ax(2),ctd(i).S_mean,ctd(i).depth_mean,'-','color',dep_clrs(i,:),'linewidth',lw);
end
% rho
for i = 1:3
    plot(ax(3),ctd(i).rho_mean,ctd(i).depth_mean,'-','color',dep_clrs(i,:),'linewidth',lw);
end
% N
for i = 1:3
    h(i) = plot(ax(4),ctd(i).N,ctd(i).depth_mean,'-','color',dep_clrs(i,:),'linewidth',lw);
end

% axes
ylabel(ax(1),'depth [m]','fontsize',fs)
for i = 1:length(ax)
    axis(ax(i),'ij')
    box(ax(i),'on')
    if i > 1
        ax(i).YTickLabel = {};
    end
    linkaxes(ax,'y')
    xlabel(ax(i),x_lbls{i},'fontsize',fs)
end
ylim(ax(1),[0 160])

xlim(ax(1),[4 8])
xlim(ax(2),[23 32])
xlim(ax(3),[1018 1025])
xlim(ax(4),[0 7e-3])

for i = 1:length(ax(4).XTickLabel)
    if str2double(ax(4).XTickLabel{i})
        ax(4).XTickLabel{i} = [ax(4).XTickLabel{i} 'x10^{-3}'];
    end
end

% text(ax(4),0.85,-0.055,'x10^{-3}','fontsize',ax(4).FontSize,'units','normalized')

% legend
lgd = legend(ax(4),h,{'21 m','44 m','52 m'},'fontsize',fs-1);
lgd.Position(1:2) = [0.815 0.4];