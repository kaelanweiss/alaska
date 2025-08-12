% Script to generate table 6 for the melt rate paper. This is a summary of
% modeled melt rate error factors for all deployments. The script generates
% everything within the \begin{tabular} ... \end{tabular} commands.
%
% Updated for Sch22 coarse 20 Mar 2025
%
% KJW
% 5 Nov 2024

clear

addpath('..')

% load data
ms_tbl = loadMSInfo(1:19,'manualwindows');
load ../../../data/melt_models.mat

% extract fields
m = msTable2Vector(ms_tbl.m);
m_ci = ms_tbl.m_ci;

% convert melt rate to m/day
m = m*.24;
m_ci = m_ci*.24;
m(end) = nan;

% error factor
r = repmat(m,[1 size(m_mean,2)])./m_mean;
r_ci = repmat(m_ci,[1 size(m_mean,2)])./m_mean;

r_mean = mean(r,'omitnan');
r_mean_ci = std(r,'omitnan');
r_med = median(r,'omitnan');

% round for printing accurately
m = round(m,2);
m_ci = round(m_ci,2);
r = [round(r(:,1:2),1) round(r(:,3),2) round(r(:,4:5),1) round(r(:,6),2)];
r_ci = [round(r_ci(:,1:2),1) round(r_ci(:,3),2) round(r_ci(:,4:5),1) round(r_ci(:,6),2)];

r_mean = [round(r_mean(1:2),1) round(r_mean(3),2) round(r_mean(4:5),1) round(r_mean(6),2)];
r_mean_ci = [round(r_mean_ci(1:2),1) round(r_mean_ci(3),2) round(r_mean_ci(4:5),1) round(r_mean_ci(6),2)];
r_med = [round(r_med(1:2),1) round(r_med(3),2) round(r_med(4:5),1) round(r_med(6),2)];

% rows
rows = 1:31;
nrows = length(rows);

% column labels and units
col_lbls = {'Segment','$m_{\\mathrm{obs}}$ [m/day]',...
    '$\\frac{m_\\mathrm{obs}}{m_\\mathrm{3eqn}}$',...
    '$\\frac{m_\\mathrm{obs}}{m_\\mathrm{3eqn,coarse}}$',...
    '$\\frac{m_\\mathrm{obs}}{m_\\mathrm{bulk}}$',...
    '$\\frac{m_\\mathrm{obs}}{m_\\mathrm{3eqn,Sch22}}$'};
col_units = {'','m/day','','','',''};

% cases
segments = [8 3 29 9];
seg_lbls = upper({'a','b','c','d'});

% header lines
ncols = length(col_lbls);
line1 = '';
line2 = '';
for i = 1:ncols
    line1 = [line1 sprintf('%s',col_lbls{i})];
    if ~isempty(col_units{i})
        line2 = [line2 sprintf('(%s)',col_units{i})];
    end
    if i<ncols
        line1 = [line1 ' & '];
        line2 = [line2 ' & '];
    end
end
line1 = ['\t\t' line1 ' \\\\\n'];
line2 = ['\t\t' line2 ' \\\\\n'];

fid = fopen('../../../../../papers/weiss_et_al_2024/table_error_latex.txt','w');
fprintf(fid,line1);
% fprintf(fid,line2);
fprintf(fid,'\t\t\\hline\n');

% create rows
for i = 1:nrows
    % data
    if ~any(i==segments)
        line = sprintf('%d & $%.2f\\pm%.2f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.2f\\pm%.2f$ \\\\\n',...
            i,m(i),m_ci(i),r(i,1),r_ci(i,1),r(i,5),r_ci(i,5),r(i,4),r_ci(i,4),r(i,6),r_ci(i,6));
    else
        line = sprintf('%d (%s) & $%.2f\\pm%.2f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.2f\\pm%.2f$ \\\\\n',...
            i,seg_lbls{i==segments},m(i),m_ci(i),r(i,1),r_ci(i,1),r(i,5),r_ci(i,5),r(i,4),r_ci(i,4),r(i,6),r_ci(i,6));
    end
    
    % don't touch this
    line = strrep(line,'pm0.','pm.');
    line = strrep(line,'\','\\');
    fprintf(fid,['\t\t' line]);
end
fprintf(fid,'\t\t\\hline\n');

% mean error factor
line = sprintf('\t\t Mean &  & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.2f\\pm%.2f$ \\\\\n',...
    r_mean(1),r_mean_ci(1),r_mean(5),r_mean_ci(5),r_mean(4),r_mean_ci(4),r_mean(6),r_mean_ci(6));
line = strrep(line,'pm0.','pm.');
line = strrep(line,'\','\\');
fprintf(fid,line);

% median error factor
line = sprintf('\t\t Median ($r$) & & $%.1f$ & $%.1f$ & $%.1f$ & $%.2f$ \\\\\n',...
    r_med(1),r_med(5),r_med(4),r_med(6));
line = strrep(line,'\','\\');
fprintf(fid,line);

fprintf(fid,'\t\t\\hline');

fclose(fid);





