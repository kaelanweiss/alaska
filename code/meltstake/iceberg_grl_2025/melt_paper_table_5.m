% Script to generate table 5 for the melt rate paper. This is a summary of
% observed and modeled melt rate for all deployments. The script generates
% everything within the \begin{tabular} ... \end{tabular} commands.
%
% Updated for Sch22 coarse 20 Mar 2025
%
% KJW
% 17 Apr 2024

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
m = round(m,2);
m_ci = round(m_ci,2);

% rows
rows = 1:31;
nrows = length(rows);

% column labels and units
col_lbls = {'Segment','$m_\\mathrm{obs}$','$m_\\mathrm{3eqn}$','$m_\\mathrm{3eqn,coarse}$','$m_\\mathrm{bulk}$','$m_\\mathrm{3eqn,Sch22}$'};
col_units = {'','m/day','m/day','m/day','m/day','m/day'};

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
        line2 = [line2 sprintf('[%s]',col_units{i})];
    end
    if i<ncols
        line1 = [line1 ' & '];
        line2 = [line2 ' & '];
    end
end
line1 = ['\t\t' line1 ' \\\\\n'];
line2 = ['\t\t' line2 ' \\\\\n'];

fid = fopen('../../../../../papers/weiss_et_al_2024/table_melt_rates_latex.txt','w');
fprintf(fid,line1);
fprintf(fid,line2);
fprintf(fid,'\t\t\\hline\n');

% create rows
for i = 1:nrows
    % data
    if ~any(i==segments)
        line = sprintf('%d & $%.2f\\pm%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\\n',...
            i,m(i),m_ci(i),m_mean(i,1),m_mean(i,5),m_mean(i,4),m_mean(i,6));
    else
        line = sprintf('%d (%s) & $%.2f\\pm%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\\n',...
        i,seg_lbls{i==segments},m(i),m_ci(i),m_mean(i,1),m_mean(i,5),m_mean(i,4),m_mean(i,6));
    end
    
    % don't touch this
    line = strrep(line,'pm0.','pm.');
    line = strrep(line,'\','\\');
    fprintf(fid,['\t\t' line]);
end
fprintf(fid,'\t\t\\hline');

fclose(fid);





