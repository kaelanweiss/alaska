% Script to plot thermistor and ADCP orientation (and eventually velocity?)
% for a given Spider deployment.
% KJW
% 29 Aug 2022
clear;
addpath('C:/Users/kweis/Desktop/School/Matlab/rbr-rsktools/');

depname = 'F:/Alaska2022/data/raw/iceberg_surveys/20220824_teaparty/spider/';

% load RBR data
rbr = parseRBRDeployment(depname);

% load Nortek data
fprintf('Nortek\n');
ntk = getNortekOrientation(depname);
ntk.rol(ntk.rol<0) = ntk.rol(ntk.rol<0)+360;

% load deployment time
tlim = [NaN NaN];
fid = fopen(fullfile(depname,'deployment_times.txt'),'r');
for i = 1:2
    line = fgetl(fid);
    while line(1) == '#'
        line = fgetl(fid);
    end
    vals = strsplit(line,',');
    tlim(i) = datenum(strjoin(vals(1:6),'-'),'yyyy-mm-dd-HH-MM-SS');
end
fclose(fid);

ntk.tidx = ntk.t>=tlim(1) & ntk.t<=tlim(2);



%% plot
t_offset = 0; % ADT is 8 hours behind UTC
ax = [];
% temperature and pressure
figure(1); clf;
hold on;

% rbr
n_rbr = length(rbr);
lbls = {};
for i = 1:n_rbr
    plot(rbr(i).time+t_offset,rbr(i).values);
    lbls{i} = num2str(rbr(i).sn);
end
legend(lbls);
ylabel('temp [^\circC]');
ax(1) = gca;

% ntk
yyaxis right;
set(gca,'YColor','k');
plot(ntk.t+t_offset,ntk.p,'k');
ylabel('pressure [dbar]');
xlim(tlim+t_offset);

grid;
datetick('x','mmmdd HH:MM','keeplimits');

% adcp orientation
figure(2); clf;
flds = {'hdg','ptc','rol'};
for i = 1:length(flds)
    ax(end+1) = subplot(length(flds),1,i);
    plot(ntk.t+t_offset,ntk.(flds{i}),'k','linewidth',1);
    grid;
    ylabel(flds{i});
    ylim(extrema(ntk.(flds{i})(ntk.tidx)));
    xlim(tlim+t_offset);
    datetick('x','mmmdd HH:MM');
end
linkaxes(ax,'x');

