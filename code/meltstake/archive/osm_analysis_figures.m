 % Script to make main results scatter plot(s) for OSM.

set(0,'defaulttextinterpreter','latex')

% load data
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');

% cm/hr to m/dy
cmh2md = 0.24;

% extract fields
z = ms_tbl.depth;
[m,m_ci] = msTable2Vector(ms_tbl.m);
[m1,m1_ci] = msTable2Vector(ms_tbl.m1);
[m2,m2_ci] = msTable2Vector(ms_tbl.m2);
[m3,m3_ci] = msTable2Vector(ms_tbl.m3);
[T,T_ci] = msTable2Vector(ms_tbl.T);
[S,S_ci] = msTable2Vector(ms_tbl.S);

% better melt error bars
m_all = [m1 m2 m3];
m_diff = abs([m1 m2 m3] - repmat(m,[1 3]));
m_diff_mean = mean(m_diff,'all','omitnan');
m_diff(isnan(m_diff)) = m_diff_mean;
m_diff_max = max(m_diff,[],2,'omitnan');
% regression uncertainty, disagreement between beams (including missing beams), ice obliquity,
% uncertainty from tank experiment
m_ci = sqrt(m_ci.^2 + m_diff_max.^2 + (m*(1-cosd(25))).^2 + (.04*m).^2);

% convert melt rate to m/day
m = m*cmh2md;
m_ci = m_ci*cmh2md;

nu = 15;
u = nan(size(m,1),nu);
u_ci = nan(size(m,1),nu);
for i = 1:nu
    [u(:,i),u_ci(:,i)] = msTable2Vector(ms_tbl.(['u' num2str(i)]));
end

% % choose u scale
% u_num = 5;
% u_plt = abs(u(:,u_num));
% %u_plt = sqrt(u(:,10).^2 + (2*u(:,11)).^2);
% %u_plt = u(:,10);
% %u_plt = 0.1*ones(size(u(:,1)));
% %u_plt = max(u(:,10),u(:,11));
% u_plt_ci = u_ci(:,u_num);

[u_plt,u_plt_ci] = msTable2Vector(ms_tbl.u0);

% standard error
tau_decorr = ms_tbl.tau_decorr;
NT = 60*ms_tbl.Duration./tau_decorr;
Nu = 60*ms_tbl.Duration./tau_decorr;

T_std = T_ci;
u_std = u_plt_ci;
S_std = S_ci;

T_ci = T_ci./sqrt(NT);
S_ci = S_ci./sqrt(4);
u_plt_ci = u_plt_ci./sqrt(Nu);

% velocity ratio
uw_ratio = abs(u(:,13))./u(:,14); % u/w
uw_ratio(u(:,14)<0) = nan;

% flow type
idx_plume = strcmp(ms_tbl.flow,'plume');
idx_ambient = strcmp(ms_tbl.flow,'ambient');
idx_wave = strcmp(ms_tbl.flow,'wave');
idx_bubble = strcmp(ms_tbl.flow,'bubble');

% flow magnitude
vel_lim = .04;
% idx_conv = u_plt<=vel_lim;
idx_conv = abs(u_plt)<=vel_lim;
% idx_conv = idx_conv & T<6;

% flow anisotropy
idx_uw = uw_ratio < 0.8;

% simple 3eqn model (S=0)
% prefactor
cw = 4.2e3;
Cd = (2.5e-3);
GammaT = (2.2e-2);
L = 3.3e5;
C = cw*sqrt(Cd)*GammaT/L; % 1/C
m_simple = C*u_plt.*T*3600*100;
m_simple = m_simple*cmh2md;

% full model
[m_3eqn,Tb,Sb] = solve3Eqn(u_plt,T,S,z);
m_3eqn = m_3eqn*(24*60*60);

% full model error bars
m_3eqn_ci = C*sqrt((T-Tb).^2.*(u_plt_ci).^2 + u_plt.^2.*(T_ci).^2);
m_3eqn_ci = m_3eqn_ci*(24*60*60);

% sort convective by T
[~,sidx] = sort(T);
sidx = intersect(sidx,find(idx_conv),'stable');

% K&M 2015 model
Td = T(sidx)-Tb(sidx);
m_KM15 = 0.25*Td.^(4/3)*(86400)/1e6; % um/s --> m/day
m_KM_error = (m(sidx)-m_KM15)./m_KM15;

% Schulz
m_Sch22 = (1e-3)*(cw/L)*Td*(24*3600); % m/s --> m/day

% Jenkins 2010 with U = 0.04;
[m_3eqn_u4,~,~] = solve3Eqn(vel_lim*ones(size(u_plt(sidx))),T(sidx),S(sidx),z(sidx));
m_3eqn_u4 = m_3eqn_u4*(24*3600);

% gammaT KW (Cb = gammaT*cw/L)
% gammaT = regress(m(sidx),(cw/L)*Td.^(4/3))/(24*3600); % m/s * C^(-1/3)
[Cb,Cb_int] = regress(m(sidx)/(24*3600),Td.^(4/3)); % m/s * C^(-4/3)
m_KW = Cb*Td.^(4/3)*(24*3600);
gammaT = Cb*L/cw;

% melt error
m_error = (m-m_3eqn)./m_3eqn;
m_error_ci = m_ci./m_3eqn;
%err_lbl = {'$\underline{m_{obs}-m_{3eqn}}$','$m_{3eqn}$'};
err_lbl = {'m$_{error}$'};

% conv melt error
m_err_conv = (m(sidx)-m_KW)./m_KW;
m_err_conv_ci = m_ci(sidx)./m_KW;

%% plot
fs = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_ex = [.2 .5 1 2];
marker_scale = 200;
ldg_lbls = cell(1,length(m_ex));
dy_lgd = 0.4;
% figure 1 parameter space (u vs T with colors corresponding to S)
figure(1); clf
set(gca,'fontsize',fs-2)
hold on
% errorbar(S,T,-T_ci,T_ci,-S_ci,S_ci,'.','color',0.5*[1 1 1])
% scatter(S,T,40,u_plt,'o','filled','markeredgecolor','k')
scatter(u_plt,T,m*marker_scale,S,'o','filled','markeredgecolor','k')
errorbar(u_plt,T,-T_std,T_std,-u_std,u_std,'linestyle','none','color',0.5*[1 1 1],'linewidth',1)

cbar = colorbar;
box on

%xlim([0 .13])
ylim([3 11])
set(gca,'clipping','off')

% homemade legend entries
clear h
XLIM = xlim;
YLIM = ylim;
for i = 1:length(m_ex)
    xi = XLIM(1)+0.75*diff(XLIM);
    yi = YLIM(1)+0.45*diff(YLIM)+(i-1)*dy_lgd;
    h(i) = scatter(xi,yi,marker_scale*m_ex(i),'ko');
    text(xi+.007,yi,sprintf('%.1f m/day',m_ex(i)),'fontsize',fs-4,'interpreter','latex')
end
text(xi+.007,yi+dy_lgd,'Melt rate','fontsize',fs-4,'interpreter','latex')

% [~,lgd2] = legend(h,lgd_lbls);
% for i = 1:length(m_ex)
%     patchi = lgd2(length(m_ex)+i).Children;
%     patchi.MarkerSize = m_ex(i)*500;
% end

xlabel('U (m/s)','fontsize',fs)
ylabel('T (C)','fontsize',fs)
cbar.Label.String = 'S (psu)';
cbar.Label.FontSize = fs;
cbar.Label.Interpreter = 'latex';
cmocean('hal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % figure 2 (m vs uT)
% figure(2); clf
% hold on
% 
% uT = u_plt.*T;
% errorbar(uT,m,m_ci,'ko','markersize',4,'markerfacecolor',colors(1))
% plot(uT,m_simple,'k^','markersize',4,'markerfacecolor','k')
% plot(uT,m_3eqn,'ko','markersize',4,'markerfacecolor','r')
% 
% box on
% 
% xlabel('u*T (m*C/s)','fontsize',fs)
% ylabel('melt rate (m/dy)','fontsize',fs)
% 
% legend({'m_{obs}','m_{3eqn} (S=0)','m_{3eqn}'},'fontsize',fs-2,'location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 3 (m vs m)
figure(3); clf
set(gca,'fontsize',fs-2)
hold on
h1 = plot([0 1.4],[0 1.4],'k-');
h2 = plot([0 1],3*[0 1],'k--');

errorbar(m_3eqn(~idx_conv),m(~idx_conv),-m_ci(~idx_conv),m_ci(~idx_conv),-m_3eqn_ci(~idx_conv),m_3eqn_ci(~idx_conv),'linestyle','none','color',0.5*[1 1 1],'linewidth',1);
errorbar(m_3eqn(idx_conv),m(idx_conv),-m_ci(idx_conv),m_ci(idx_conv),-m_3eqn_ci(idx_conv),m_3eqn_ci(idx_conv),'linestyle','none','color',0.5*[1 1 1],'linewidth',1);

h3 = plot(m_3eqn(idx_conv),m(idx_conv),'o','markersize',5,...
    'markerfacecolor',colors(1),'markeredgecolor','k');
h4 = plot(m_3eqn(~idx_conv),m(~idx_conv),'^','markersize',5,...
    'markerfacecolor',colors(2),'markeredgecolor','k');

box on
xlim([0 max(ceil(5*m_3eqn))/5])
ylim([0 2.5])

xlabel('m (model) (m/day)','fontsize',fs)
ylabel('m (observed) (m/day)','fontsize',fs)

legend([h3 h4],{sprintf('U$\\le$%.2f m/s',vel_lim),sprintf('U$>$%.2f m/s',vel_lim)},...
    'fontsize',fs-2,'location','southeast','interpreter','latex')
% legend({'1:1'},'fontsize',fs-2,'location','southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 4 (m-m vs u,T,S)
figure(4); clf
x_data = {u_plt,T,S};
x_ci = {u_plt_ci,T_ci,S_ci};
x_lbls = {'U (m/s)','T (C)','S (psu)'};
clrs = [colors(1); colors(3); colors(2)];
for i = 1:3
    axes(figure(4),'position',axgridpos(3,1,i,.08,.08,.12,.1),'fontsize',fs-2)
    hold on
    clear h
    
    if i==1
%         plot(sort(x_data{i}),-1+(gammaT/(sqrt(CD)*0.5*GammaT)).*sort(x_data{i}).^-1,'k--')
        plot(sort(x_data{i}),-1+((Cb*L/cw)/(sqrt(Cd)*0.5*GammaT)).*sort(x_data{i}).^-1,'k--')
        plot(sort(x_data{i}),-1+((Cb_int(1)*L/cw)/(sqrt(Cd)*0.5*GammaT)).*sort(x_data{i}).^-1,'k--')
        plot(sort(x_data{i}),-1+((Cb_int(2)*L/cw)/(sqrt(Cd)*0.5*GammaT)).*sort(x_data{i}).^-1,'k--')
%         plot(sort(x_data{i}),-1+(((4.4e-7)*L/cw)/(sqrt(Cd)*0.5*GammaT)).*sort(x_data{i}).^-1,'r--')
%         plot(sort(x_data{i}),-1+(((7.8e-7)*L/cw)/(sqrt(Cd)*0.5*GammaT)).*sort(x_data{i}).^-1,'r--')
        plot(sort(x_data{i}),-1+(((2.5e-7)*L/cw)/(sqrt(Cd)*0.5*GammaT)).*sort(x_data{i}).^-1,'g--')
        plot(0.04*[1 1], [-1 7],'k:')
        plot([0 .12],0*[1 1],'k-')
    end

    % 3eqn
    %errorbar(x_data{i},(m)./m_3eqn,m_ci,'ko','markersize',4,'markerfacecolor',clrs(i,:))
    %errorbar(x_data{i},m_error,-m_error_ci,m_error_ci,-x_ci{i},x_ci{i},'linestyle','none','color',0.5*[1 1 1],'linewidth',1)
    errorbar(x_data{i},m_error,m_error_ci,'linestyle','none','color',0.5*[1 1 1],'linewidth',1)
    h(1) = plot(x_data{i},m_error,'ko','markerfacecolor',colors(i),'markersize',4.5);
    
    %conv
    errorbar(x_data{i}(sidx),m_err_conv,m_err_conv_ci,'linestyle','none','color',0.5*[1 1 1],'linewidth',1)
    h(2) = plot(x_data{i}(sidx),m_err_conv,'ko','markerfacecolor','none','markersize',4.5);
    
    xlabel(x_lbls{i},'fontsize',fs)
    ylabel(err_lbl,'fontsize',fs)
    box on

    if i==1
        legend(h,{'shear','convective'},'orientation','horizontal','location','northeast','fontsize',fs-2)
    end

    ylim([-1 ceil(max(m_error))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 5 (m vs u,T,S)
figure(5); clf
x_data = {u_plt,T,S};
x_ci = {u_plt_ci,T_ci,S_ci};
x_lbls = {'U (m/s)','T (C)','S (psu)'};
clrs = [colors(1); colors(3); colors(2)];
for i = 1:3
    axes(figure(5),'position',axgridpos(3,1,i,.08,.08,.12,.1))
    hold on
    errorbar(x_data{i},m,-m_ci,m_ci,-x_ci{i},x_ci{i},'linestyle','none','color',0.5*[1 1 1],'linewidth',1)
    plot(x_data{i},m,'ko','markersize',5,'markerfacecolor',colors(i))
    plot(x_data{i}(sidx),m(sidx),'ko','markersize',5,'markerfacecolor','k')
    xlabel(x_lbls{i},'fontsize',fs)
    ylabel('m (m/day)','fontsize',fs)
    ylim([0 2.5])
    box on
    set(gca,'fontsize',fs-2)
    if i==3
        xlim([22 29])
        set(gca,'clipping','off')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % figure 5 (m error vs w/u ratio)
% figure(5); clf
% plot(uw_ratio,m_error,'ko','markerfacecolor',colors(1),'markersize',4)
% 
% ylabel(err_lbl,'fontsize',fs)
% xlabel('u/w','fontsize',fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 6 (m vs m with flow type distinction)
figure(6); clf
hold on
plot([0 2],[0 2],'k-')
errorbar(m_3eqn,m,m_ci,'.','markersize',4,'color',0.5*[1 1 1],'linewidth',1)

plot(m_3eqn(idx_ambient),m(idx_ambient),'ko','markersize',6,'markerfacecolor',colors(1))
plot(m_3eqn(idx_plume),m(idx_plume),'k^','markersize',6,'markerfacecolor',colors(2))
plot(m_3eqn(idx_wave),m(idx_wave),'ksquare','markersize',6,'markerfacecolor',colors(3))
plot(m_3eqn(idx_bubble),m(idx_bubble),'khexagram','markersize',8,'markerfacecolor',colors(4))

box on
xlim([0 2])
ylim([0 2])

xlabel('predicted melt rate (m/day)','fontsize',fs)
ylabel('observed melt rate (m/day)','fontsize',fs)

legend({'1:1','CI','ambient','plume','wave','bubble'},'fontsize',fs-2,'location','southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % figure 7 (m error vs uTS with flow type distinction)
% figure(7); clf
% x_data = {u_plt,T,S};
% x_lbls = {'u (m/s)','T (C)','S (psu)'};
% clrs = [colors(1); colors(3); colors(2)];
% for i = 1:3
%     axes(figure(7),'position',axgridpos(3,1,i,.08,.08,.12,.1))
%     hold on
%     plot(x_data{i}(idx_ambient),m_error(idx_ambient),'ko','markersize',6,'markerfacecolor',colors(1))
%     plot(x_data{i}(idx_plume),m_error(idx_plume),'k^','markersize',6,'markerfacecolor',colors(2))
%     plot(x_data{i}(idx_wave),m_error(idx_wave),'ksquare','markersize',6,'markerfacecolor',colors(3))
%     plot(x_data{i}(idx_bubble),m_error(idx_bubble),'khexagram','markersize',8,'markerfacecolor',colors(4))
% 
%     xlabel(x_lbls{i},'fontsize',fs)
%     ylabel(err_lbl,'fontsize',fs)
%     if i == 1
%         legend({'ambient','plume','wave','bubble'},'fontsize',fs-2,'location','best')
%     end
% 
%     box on
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % figure 8 (m vs m with flow > 0.04 distinction)
% figure(8); clf
% hold on
% plot([0 1.2],[0 1.2],'k-')
% errorbar(m_3eqn,m,m_ci,'.','markersize',4,'color',0.5*[1 1 1],'linewidth',1)
% 
% plot(m_3eqn(~idx_conv),m(~idx_conv),'ko','markersize',6,'markerfacecolor',colors(1))
% plot(m_3eqn(idx_conv),m(idx_conv),'k^','markersize',6,'markerfacecolor',colors(2))
% 
% 
% box on
% ylim([0 2.5])
% 
% xlabel('predicted melt rate (m/day)','fontsize',fs)
% ylabel('observed melt rate (m/day)','fontsize',fs)
% 
% legend({'1:1','CI',sprintf('U>%.1f cm/s',100*vel_lim),sprintf('U<=%.1f cm/s',100*vel_lim)},'fontsize',fs-2,'location','southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % figure 9 (m vs m with flow anisotropy distinction)
% figure(9); clf
% hold on
% plot([0 1.2],[0 1.2],'k-')
% errorbar(m_3eqn,m,m_ci,'.','markersize',4,'color',0.5*[1 1 1],'linewidth',1)
% 
% plot(m_3eqn(~idx_uw),m(~idx_uw),'ko','markersize',6,'markerfacecolor',colors(1))
% plot(m_3eqn(idx_uw),m(idx_uw),'k^','markersize',6,'markerfacecolor',colors(2))
% 
% box on
% ylim([0 2.5])
% 
% xlabel('predicted melt rate (m/day)','fontsize',fs)
% ylabel('observed melt rate (m/day)','fontsize',fs)
% 
% legend({'1:1','CI','w<u','w>u'},'fontsize',fs-2,'location','southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 10 (m vs T in convective regime)
figure(10); clf
hold on

h1 = plot(Td,m_KW,'k--','linewidth',1);
h2 = plot(Td,m_3eqn_u4,'k-','linewidth',1);
h3 = plot(Td,m_KM15,'-','linewidth',1,'color',colors(5));
h4 = plot(Td,m_Sch22,'-','linewidth',1,'color',colors(7));

errorbar(T(sidx)-Tb(sidx),m(sidx),-m_ci(sidx),m_ci(sidx),-T_ci(sidx),T_ci(sidx),'linestyle','none','color',0.5*[1 1 1],'linewidth',1)

h5 = plot(Td,m(sidx),'ko','markerfacecolor',colors(1),'markersize',5);

xlabel('$T-T_b$ (C)','fontsize',fs)
ylabel('m (m/day)','fontsize',fs)

box on
set(gca,'yscale','log')
grid on
ylim([1e-1, 1e1])

legend([h1 h5 h2 h3 h4],{sprintf('fit (C=%.1e)',Cb),'obs (U$\le$.04m/s)','3eqn (U=.04m/s)','K\&M 2015, Gayen 2016','Schulz 2022'},...
    'location','northoutside','interpreter','latex','fontsize',fs-2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % figure 11 (m-m vs u,T,S in convective regime)
% figure(11); clf
% x_data = {u_plt,T,S};
% x_lbls = {'u (m/s)','T (C)','S (psu)'};
% clrs = [colors(1); colors(3); colors(2)];
% for i = 1:3
%     axes(figure(11),'position',axgridpos(3,1,i,.08,.08,.12,.1))
%     %errorbar(x_data{i},(m)./m_3eqn,m_ci,'ko','markersize',4,'markerfacecolor',clrs(i,:))
%     plot(x_data{i}(sidx),m_KM_error,'ko','markersize',4,'markerfacecolor',clrs(i,:))
%     xlabel(x_lbls{i},'fontsize',fs)
%     ylabel(err_lbl,'fontsize',fs)
%     ylim([0 8])
% end
