% Script to calculate melt rates using 881a sonar data.
%
% KJW
% 4 Feb 2026
clear

proc_dir = 'F:/meltstake/data/proc';

seg_tbl = loadMSInfo(26:28,'segments');
dep_names = unique(seg_tbl.Folder,'stable');
ndeps = length(dep_names);

% load
for i = ndeps:-1:1
    load_edges = load(fullfile(proc_dir,dep_names{i},'sonar881a_edges.mat'));
    edges(i) = load_edges.sonar881a_edges;
end

% QUICK FIX FOR LAGGING 881 TIME AXIS
edges(2).time = edges(2).time + minutes(1);

%% melt rates
% recipe (assuming no ice angle)
%   1. split edges into segments
%   2. median/sigma filter
%   3. linear regression at each angle
%   4. combine regressions with var-wgtd mean
m = nan(size(seg_tbl,1),1);
m_ci = nan(size(m));
pos = 1;
for i = 2
    segments = seg_tbl(strcmp(seg_tbl.Folder,dep_names{i}),:);
    nsegs = size(segments,1);
    for j = 7%1:nsegs
        % time index
        t1 = segments.Start(j);
        t2 = segments.End(j);
        idxt = edges(i).time>=t1 & edges(i).time<=t2;
        tj = edges(i).time(idxt);
        % angle index
        idxa = abs(edges(i).angle)<=30;
        x = edges(i).x_edge(idxt,idxa);
        r = edges(i).r_edge(idxt,idxa);
        % clean up using standard deviation
        std_lim = 1.5*1*days(tj(end)-tj(1));
        idx_std = std(x,0,1) <= std_lim;
        x = x(:,idx_std);
        % regression
        [nt,na] = size(x);
        fit_params = nan(na,2,3); % two parameters (m,b) and 3 vals (actual and conf intervals)
        for k = 1:na
            if nt < 2
                break
            end
            [b,bint] = regress(x(:,k),[ones(nt,1) days(tj-tj(1))]);
            if all(isnan(bint),'all')
                bint(:) = 1;
            end
            fit_params(k,:,:) = [flip(b) flip(bint)];
        end
        mj = fit_params(:,1,1);
        dmj = diff(squeeze(fit_params(:,1,2:3)),[],2);
        m(pos) = varWgtMean(mj,dmj);
        fprintf('%d.%d: %.2f %.3f\n',i,j,round(m(pos),2),round(na/length(idxa),3))
        pos = pos+1;
        % plot for debugging
        figure(100*i+10*j); clf
        hold on
        plot(tj,x,'.-')
    end
    fprintf('\n')
end




