% Script to plot a bunch of data records to line up time axes
% KJW
% 30 Aug 2022

clear;

raw_path = 'F:/Alaska2022/data/iceberg_surveys/raw/';
berg = '20220825_singingflower';
deppath = fullfile(raw_path,berg,'spider');

%% ADCP
fprintf('ADCP...\n');
ntk = getNortekOrientation(deppath);

%% ROV
% fprintf('ROV...\n');
% rov = parseROVTelemetry(deppath);

%% RBR (could turn this into its own function)
fprintf('RBR...\n');
% read instruments file
rbr = parseRBRDeployment(deppath);
n_rbr = length(rbr);

%% Temperature from solos and pressure from ADCP
ax = [];
figure(1); clf;
lbls1 = cell(n_rbr,1);
subplot(6,1,3:6); hold on;
ax(end+1) = gca;
for i = 1:n_rbr
    plot(rbr(i).time,rbr(i).values(:,1),'.-','linewidth',1.2);
    lbls1{i} = num2str(rbr(i).sn);
end
%xlim([datenum(2022,8,25,17,13,20) datenum(2022,8,25,17,13,25)]);
grid;
datetick('x','mmmdd HH:MM:SS','keeplimits');
lgd = legend(lbls1,'fontsize',12);
ylabel('temperature');
set(gca,'clipping','off');

subplot(6,1,1:2);
ax(end+1) = gca;
plot(ntk.t-0/86400,ntk.ptc,'k');
grid;
ylabel('pressure');
set(gca,'clipping','off');

linkaxes(ax,'x');
datetick('x','mmmdd HH:MM:SS','keeplimits');

%% Pressure from ADCP, ROV, concerto, and duet(s)
% figure(2); clf; hold on;
% lbls2 = {'duet','concerto','adcp','rov'};
% % duet
% plot(rbr(13).time,rbr(13).values(:,2),'k.-','linewidth',1.2);
% 
% % concerto
% plot(rbr(14).time,rbr(14).values(:,3)+0.22,'.-','linewidth',1.2);
% 
% % adcp
% plot(ntk.t-(6/86400),ntk.p+10.025,'.-','linewidth',1,'Color',colors(4));
% 
% % rov
% plot(rov.time-(3/86400),-rov.altitudeAMSL+11.01,'.-','linewidth',1,'Color',colors(5));
% 
% xlim([datenum(2022,8,25,17,15,0) datenum(2022,8,25,19,0,0)]);
% grid;
% datetick('x','mmmdd HH:MM:SS.fff','keeplimits');
% legend(lbls2,'fontsize',12);

