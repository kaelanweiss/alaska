% Script to plot a bunch of data records to line up time axes
% KJW
% 30 Aug 2022

clear;

deppath = 'F:/Alaska2022/data/iceberg_surveys/raw/20220821_marchhare/dragon';

%% ADCP
fprintf('ADCP...\n');
ntk = getNortekOrientation(deppath);

%% ROV
fprintf('ROV...\n');
rov = parseROVTelemetry(deppath);

%% RBR (could turn this into its own function)
fprintf('RBR...\n');
% read instruments file
rbr = parseRBRDeployment(deppath);

%% Temperature from solos and duet(s)
figure(1); clf; hold on;
lbls1 = cell(length(rbr),1);
for i = 1:7
    plot(rbr(i).time,rbr(i).values(:,1),'.-','linewidth',1.2);
    lbls1{i} = num2str(rbr(i).sn);
end
for i = 8:length(rbr)
    plot(rbr(i).time,rbr(i).values(:,1),'.--','linewidth',1.2);
    lbls1{i} = num2str(rbr(i).sn);
end
%xlim([datenum(2022,8,25,17,13,20) datenum(2022,8,25,17,13,25)]);
grid;
datetick('x','mmmdd HH:MM:SS.fff','keeplimits');
lgd = legend(lbls1,'fontsize',12);
title('Rake temperature');

%% Pressure from ADCP, ROV, concerto, and duet(s)
figure(2); clf; hold on;
%lbls2 = {'duet','concerto','adcp','rov'};
lbls2 = {};
% duet
duet_idx = find(strcmp('duet',{rbr.model}),1);
if ~isempty(duet_idx)
    plot(rbr(duet_idx).time,rbr(duet_idx).values(:,2),'k.-','linewidth',1.2);
    lbls2{end+1} = 'duet';
end

% concerto
concerto_idx = find(strcmp('concerto',{rbr.model}),1);
if ~isempty(concerto_idx)
    plot(rbr(concerto_idx).time,rbr(concerto_idx).values(:,3)+0.22,'.-','linewidth',1.2);
    lbls2{end+1} = 'concerto';
end

% adcp
if ~isempty(fieldnames(ntk))
    plot(ntk.t-(7.25/86400),ntk.p+10.025,'.-','linewidth',1,'Color',colors(4));
    lbls2{end+1} = 'adcp';
end

% rov
if ~isempty(fieldnames(rov))
    plot(rov.time-(1/86400),medianFilter(-rov.altitudeAMSL,15)+11.01,'.-','linewidth',1,'Color',colors(5));
    lbls2{end+1} = 'rov';
end

%xlim([datenum(2022,8,25,17,15,0) datenum(2022,8,25,19,0,0)]);
grid;
datetick('x','mmmdd HH:MM:SS.fff','keeplimits');
legend(lbls2,'fontsize',12);

