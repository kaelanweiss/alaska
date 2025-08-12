% Script to calculate melt rates during meltstake deployments.
%
% KJW
% 6 Feb 2024

clear

ms_tbl = loadMSInfo('manualwindows');
ms_filt = loadMSInfo('melt_filter');

dep_nums = unique(ms_tbl.Number,'stable');
dep_names = unique(ms_tbl.Folder,'stable');

% select deployment index
dep_num = 27;

% corresponding rows and deployment name
idx_dep = dep_num==ms_tbl.Number;
dep_name = ms_tbl.Folder{find(idx_dep,1)};

% time window
idx_dep = dep_num==ms_tbl.Number;
fprintf('\n====== DEPLOYMENT %d (%s) ======\n',dep_num,dep_name)
t_start = ms_tbl.Start(idx_dep);
t_end = ms_tbl.End(idx_dep);

% load data
load(fullfile('F:meltstake\data\proc',dep_name,'adv','pck_edges.mat'))

% plot all
figure(1); clf
plot(edges.time,edges.pos,'.-')

%% window
% select window
j = 1;

fprintf('--- Window %d ---\n',j)

% find row in ms_tbl
row_num = find(ms_tbl.Number==dep_num & ms_tbl.Window==j);
if isempty(row_num)
    warning('Window %d not found',j)
    return
end

% time limits
t1 = ms_tbl.Start(row_num);
t2 = ms_tbl.End(row_num);
idxt = edges.time>=t1 & edges.time<=t2;

% adv beams
beams = str2double(strsplit(ms_tbl.Beam{row_num},','));
nb = length(beams);

% filter parameters
nsigma = [ms_filt(row_num,:).n_sigma1 ms_filt(row_num,:).n_sigma2 ms_filt(row_num,:).n_sigma3];
npass = [ms_filt(row_num,:).n_pass1 ms_filt(row_num,:).n_pass2 ms_filt(row_num,:).n_pass3];

% extract data
time = edges.time(idxt);
pos = edges.pos(idxt,beams);

figure(2); clf; hold on
for k = 1:nb
    plot(time,pos(:,k),'.-','color',colors(beams(k)))
end

% fix jumps
idxd = tdrill>=t1 & tdrill<=t2;
tdj = tdrill(idxd);

jump = [38 0 44]+1.5;
% jump = [0 26.3 0];

for k = 1:length(tdj)
    idxp = time>tdj;
    for p = 1:nb
        pos(idxp,p) = pos(idxp,p) + jump(k,beams(p));
    end
end

figure(3); clf; hold on
for k = 1:nb
    plot(time,pos(:,k),'.-','color',colors(beams(k)))
end

% remove outliers
for k = 1:nb
    [~,idxs] = sigmaFilter(detrend(pos(:,k)),nsigma(beams(k)),1,npass(beams(k)));
    plot(time(idxs),pos(idxs,k),'kx','markersize',8,'linewidth',1)
    pos(idxs,k) = nan;
end

figure(4); clf; hold on
for k = 1:nb
    plot(time,pos(:,k),'.-','color',colors(beams(k)))
end

% linear regression
tfit = hours(time-time(1));
for k = 1:nb
    [b,bint,r,rint,stats] = regress(pos(:,k)/10,[ones(size(tfit)) tfit],0.05);
    fprintf('m%d=%.1f (%.1f)\tnsigma=%.1f npass=%d\n',beams(k),round(b(2),1),round(diff(bint(2,:))/2,1),nsigma(beams(k)),npass(beams(k)))
    plot(time,10*(b(2)*tfit+b(1)),'k-')
end


