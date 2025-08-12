% Script to create figure 2 for the melt rate GRL paper. This figure has a
% panel showing the full observational parameter space as a scatter plot
% with dimensions (u, T, m, S) where the last two correspond to size and
% color. An additional 3 panels show melt rate as a function of u, T, and S
% separately.
%
% KJW
% 7 Apr 2024
set(0,'defaulttextinterpreter','latex');

clear

% load data
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

% cm/hr to m/dy
cmh2md = 0.24;

% extract fields
[m,m_ci] = msTable2Vector(ms_tbl.m);
[T,T_ci] = msTable2Vector(ms_tbl.T);
[S,S_ci] = msTable2Vector(ms_tbl.S);
[u,u_ci] = msTable2Vector(ms_tbl.u0);
u(end) = nan;

% convert melt rate to m/day
m = m*cmh2md;
m_ci = m_ci*cmh2md;

% standard error
tau_decorr = ms_tbl.tau_decorr;
NT = 60*ms_tbl.Duration./tau_decorr;
Nu = 60*ms_tbl.Duration./tau_decorr;

T_std = T_ci;
u_std = u_ci;
S_std = S_ci;

T_ci = T_ci./sqrt(NT);
S_ci = S_ci./sqrt(4);
u_ci = u_ci./sqrt(Nu);

%% plot
lw = 0.8;
fs = 11;
m_ex = [.2 .5 1 2];
marker_scale = 5;
ldg_lbls = cell(1,length(m_ex));
dy_lgd = 0.4;

pts = [8 9 15 29];
pt_lbls = upper({'c','a','b','d'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2.1 parameter space (u vs T with colors corresponding to S)
figure(21); clf
% create and size figure
figsize = [10 10];
set(figure(21),'units','centimeters'); 
fpos = get(figure(21),'position');
fpos(3:4) = figsize;
set(figure(21),'position',fpos);
set(figure(21),'paperunits','centimeters')
set(figure(21),'papersize',figsize)

% plot
ax1 = axes(figure(21));
hold on
scatter(u,T,S*marker_scale,m,'o','filled','markeredgecolor','k')
% errorbar(u,T,-T_std,T_std,-u_std,u_std,'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
text(.01,.98,'(a)','units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left')
for i = 1:length(pts)
    text(ax1,u(pts(i)),T(pts(i))+.1,sprintf('%s',pt_lbls{i}),'color','k',...
        'horizontalalignment','center','verticalalignment','bottom')
end

cbar = colorbar;
cbar.Location = 'northoutside';
box on

xlim([0 .15])
ylim([2.5 11])
clim([.24 1.8])

% % homemade legend entries
% clear h
% XLIM = xlim;
% YLIM = ylim;
% for i = 1:length(m_ex)
%     xi = XLIM(1)+0.8*diff(XLIM);
%     yi = YLIM(1)+0.45*diff(YLIM)+(i-1)*dy_lgd;
%     h(i) = scatter(xi,yi,marker_scale*m_ex(i),'ko');
%     text(xi+.015,yi,sprintf('%.1f',m_ex(i)),'interpreter','latex')
% end
% text(xi-.007,yi+2.5*dy_lgd,{'melt rate','(m/day)'},'interpreter','latex')

% labels
xlabel('U (m/s)','fontsize',fs)
ylabel('T (C)','fontsize',fs)
cbar.Label.String = 'melt rate (m/day)';
cbar.Label.FontSize = fs;
cbar.Label.Interpreter = 'latex';
colormap(flipud('hot'))

% tighten up colorbar spacing
axpos = ax1.Position;
ax1.Position = axpos + [0 0 0 .05];
cbar.Position(4) = .04;
cbar.Position(2) = sum(ax1.Position([2 4])) + .012;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 2.2 (m vs u,T,S)
figure(22); clf
% create and size figure
figsize = [9 10];
set(figure(22),'units','centimeters'); 
fpos = get(figure(22),'position');
fpos(3:4) = figsize;
set(figure(22),'position',fpos);
set(figure(22),'paperunits','centimeters')
set(figure(22),'papersize',figsize)


x_data = {u,T,S};
x_ci = {u_ci,T_ci,S_ci};
x_lbls = {'U (m/s)','T (C)','S (psu)'};
clrs = [colors(1); colors(3); colors(2)];
for i = 1:3
    axes(figure(22),'position',axgridpos(3,1,i,.1,.08,.15,.1))
    hold on
    errorbar(x_data{i},m,-m_ci,m_ci,-x_ci{i},x_ci{i},'linestyle','none','color',0.5*[1 1 1],'linewidth',lw)
    plot(x_data{i},m,'ko','markersize',4,'markerfacecolor',colors(i))
    xlabel(x_lbls{i},'fontsize',fs)
    ylabel({'melt rate','(m/day)'},'fontsize',fs)
    ylim([0 2.5])
    box on
    if i==3
        xlim([22 29])
        set(gca,'clipping','off')
    end
    text(.01,.98,['(' char(97+i) ')'],'units','normalized','fontsize',fs,...
        'verticalalignment','top','horizontalalignment','left')
end