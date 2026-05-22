clear

% melt rate from calculated distance (cm/hr)
m1 = [nan, 1.3, nan, 1.4,...
    4.7, 1.9, 5.2, nan, nan, 2.4, nan, 1.4, 5.7, nan, 2.3, 3.3];
% melt rate from beam distances (cm/hr)
m2 = [0.24, 1.0, 0.9, 1.0,...
    3.4, 1.3, 3.1, 1.5, nan, 1.5, nan, 1.3, 3.2, 1.3, 1.4, 2.4];
% temperature scale (C)
T = [nan, 3.8, 3.6 3.3,...
    9.5, 6.5, 11, 6.5, 9, 6.5, 7.5, 6, 10, 7, 8, 10];
% vel scale (m/s)
u = [.041, .019, .016, .054,...
    .046, .024, .07, .025, .058, .047, .036, .037, .049, .023, .011, .027];

% prefactor
cw = 4.2e3;
CD = (2.5e-3);
GammaT = (2.2e-2);
L = 3.3e5;
C = cw*sqrt(CD)*GammaT/L; % 1/C

% predicted melt rate
m_3eqn = C*u.*T*(100*3600);

% transfer coefficient between calculated distance melt and beam distance
% melt rate
[f,g] = fit(m1(~isnan(m1))',m2(~isnan(m1))','poly1');
b2d = 1/f.p1;
m2 = b2d*m2;

% m1 vs m2
figure(1); clf
plot(m1,m2,'ko','markerfacecolor','k')
xlabel('melt rate (calc) (cm/hr)')
ylabel('melt rate (beam dist) (cm/hr)')
box on
grid on
ylim([0 max(ylim)])
xlim([0 max(xlim)])

% T vs m2
figure(2); clf
plot(T,m2,'ko','markerfacecolor','k')
ylabel('melt rate (beam dist) (cm/hr)')
xlabel('temperature scale (C)')
box on
grid on
ylim([0 max(ylim)])
xlim([0 max(xlim)])

% u vs m2
figure(3); clf
plot(u,m2,'ko','markerfacecolor','k')
ylabel('melt rate (beam dist) (cm/hr)')
xlabel('velocity scale (m/s)')
box on
grid on
ylim([0 max(ylim)])
xlim([0 max(xlim)])

% T vs m2/u
figure(4); clf
plot(T,m2./u,'ko','markerfacecolor','k')
ylabel('m/u')
xlabel('temperature scale (C)')
box on
grid on
ylim([0 max(ylim)])
xlim([0 max(xlim)])

% T vs m2/u
figure(5); clf
plot(u,m2./T,'ko','markerfacecolor','k')
ylabel('m/T')
xlabel('velocity scale (m/s)')
box on
grid on
ylim([0 max(ylim)])
xlim([0 max(xlim)])

%%

c2m = 24/100;

figure(6); clf
% T vs m2/u
% figure(6); clf; hold on
axes(figure(6),'position',axgridpos(1,2,2,.08,.12));
hold on
plot(u.*T,m_3eqn*c2m,'color',colors(2),'linewidth',1)
plot(u(1:4).*T(1:4),m2(1:4)*c2m,'ko','markerfacecolor',0.7*[1 1 1])
plot(u(5:end).*T(5:end),m2(5:end)*c2m,'ko','markerfacecolor','k')
ylabel('melt rate (m/dy)')
xlabel('u*T (m*C/s)')
box on
grid on
ylim([0 max(ylim)])
xlim([0 max(xlim)])
legend({'melt parameterization','May','July'},'location','southeast')
text(0.02,0.96,'b)','units','normalized','fontsize',12)

axes(figure(6),'position',axgridpos(1,2,1,.08,.12));
hold on
plot(u,T,'ko','markerfacecolor','k')
plot(u(1:4),T(1:4),'ko','markerfacecolor',0.7*[1 1 1])
ylabel('T (C)')
xlabel('u (m/s)')
box on
grid on
ylim([0 12])
xlim([0 .08])
text(0.02,0.96,'a)','units','normalized','fontsize',12)

%%
figure(7); clf; hold on
plot(u.*T,m_3eqn*c2m,'k','linewidth',1)
scatter(u.*T,m2*c2m,u*7000,T,'filled','markeredgecolor','k')
ylabel('melt rate (m/dy)')
xlabel('u*T (m*C/s)')
box on
grid on
ylim([0 max(ylim)])
xlim([0 max(xlim)])
cbar = colorbar;
cbar.Label.String = 'T (C)';
cmocean('thermal')




