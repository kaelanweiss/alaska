% Script to manually set profile times given a pressure record. The profile
% times are saved as a text output.
%
% Kaelan Weiss
% 8 Sep 2022

raw_path = 'F:/Alaska2022/data/iceberg_surveys/raw/';
berg = '20220821_marchhare';

% load and find a concerto
rbr = parseRBRDeployment(fullfile(raw_path,berg,'dragon'));
concerto_idx = find(strcmp('concerto',{rbr.model}),1);
ctd = struct('time',rbr(concerto_idx).time);
ctd.P = rbr(concerto_idx).values(:,strcmp('Pressure',rbr(concerto_idx).channels));

% instantiate time vectors
t1 = []; % beginning of profiles
t2 = []; % end of profiles

% plot
figure(1); clf; hold on;
plot(ctd.time,ctd.P,'k.-');
grid;
axis ij;
datetick('x','HH:MM','keeplimits');

%% loop over this as many times as needed, changing plot limits each time
n = 4; % number of profiles
[tpts,~] = ginput(2*n);
t1 = [t1;tpts(1:2:end)];
t2 = [t2;tpts(2:2:end)];

%% save files
t1 = floor(86400*t1)/86400;
t2 = ceil(86400*t2)/86400;

% mat
save(fullfile(raw_path,berg,'dragon','profile_times.mat'),'t1','t2');

% text
fmt = 'yyyy-mm-dd HH:MM:SS.fff';
fid = fopen(fullfile(raw_path,berg,'dragon','profile_times.txt'),'w');
fprintf(fid,'#dragon profile times for %s\n',berg);
fprintf(fid,'#UTC\n');
for i = 1:length(t1)
    line = strjoin({datestr(t1(i),fmt),datestr(t2(i),fmt)},',');
    fprintf(fid,[line '\n']);
    t_check = datenum(strsplit(line,','));
    plot(t_check(1),interp1(ctd.time,ctd.P,t_check(1)),'go','markerfacecolor','g')
    plot(t_check(2),interp1(ctd.time,ctd.P,t_check(2)),'ro','markerfacecolor','r')
    text(t_check(1),interp1(ctd.time,ctd.P,t_check(1))-1,num2str(i));
    text(t_check(2),interp1(ctd.time,ctd.P,t_check(2))+1,num2str(i));
end
fclose(fid);

% plot(t1,interp1(ctd.time,ctd.P,t1),'ro','markerfacecolor','r')
% plot(t2,interp1(ctd.time,ctd.P,t2),'bo','markerfacecolor','b')