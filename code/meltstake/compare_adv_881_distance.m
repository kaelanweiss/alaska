% Script to compare ice locations in ADV and 881a sonar, mainly as a sanity
% check on ADV locations.
%
% KJW
% 22 May 2026
clear

proc_dir = 'F:/meltstake/data/proc';

% deployment number
dep_num = 28;
ms_tbl = loadMSInfo(dep_num,'segments');
nsegs = size(ms_tbl,1);
dep_name = ms_tbl.Folder{1};

% load data
load(fullfile(proc_dir,dep_name,'adv','pck_edges.mat'))
load(fullfile(proc_dir,dep_name,'sonar881a_edges.mat'))
sonar = sonar881a_edges;

%% adv location info (in sonar coordinates
y_adv = 0.1;%0.276+.00; % m
x_adv = 0.63-0.15;%-0.12; % m

% get sonar edge location corresponding to ADV
dy = .1; % m
nt = length(sonar.time);
x_881 = cell(nt,1);
for i = 1:nt
    idxi = abs(sonar.y_edge(i,:) - y_adv) <= dy;
    x_881{i} = sonar.x_edge(i,idxi);
end

%% plot
figure(1); clf
ax = axes(figure(1));
hold on

plot(edges.time,edges.pos(:,:)/1000+x_adv,'.','color',0.5*[1 1 1])
plot(edges.time,edges.boundary_location/1000+x_adv,'r-')
for i = 1:nt
    ti = repmat(sonar.time(i),size(x_881{i}));
    if length(ti) > 4
        errorbar(ti(1),mean(x_881{i}),std(x_881{i}),'k^','markersize',4,'markerfacecolor',colors(6),'linewidth',1.2)
    else
        scatter(ti,x_881{i},12,'k^','markerfacecolor',colors(6),'linewidth',1)
    end
end

ylim([0.3 0.9])
