% Script to generate table 4 for the melt rate paper. This is a summary of
% deployment data for all deployments. The script generates
% everything within the \begin{tabular} ... \end{tabular} commands.
%
% Updated to include U_coarse 20 Mar 2025
%
% KJW
% 17 Apr 2024

clear

addpath('..')

% load data
ms_tbl = loadMSInfo(1:19,'manualwindows');

% extract fields
m = msTable2Vector(ms_tbl.m);
m_ci = ms_tbl.m_ci;
[T,T_ci] = msTable2Vector(ms_tbl.T);
[S,S_ci] = msTable2Vector(ms_tbl.S);
[u,u_ci] = msTable2Vector(ms_tbl.u_max);
[uc,uc_ci] = msTable2Vector(ms_tbl.u_coarse);

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
col_lbls = {'Segment','$\\overline{U_\\mathrm{max}}$','$\\overline{U_\\mathrm{coarse}}$','$\\overline{T}$','$\\overline{S_w}$','$m_\\mathrm{obs}$'};
col_units = {'','m/s','m/s','K','psu','m/day'};

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

fid = fopen('../../../../../papers/weiss_et_al_2024/table_mUTS_latex.txt','w');
fprintf(fid,line1);
fprintf(fid,line2);
fprintf(fid,'\t\t\\hline\n');

% create rows
for i = 1:nrows
    if ~any(i==segments)
        line = sprintf('%d & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.2f\\pm%.2f$ \\\\\n',...
            i,u(i),u_ci(i),uc(i),uc_ci(i),T(i),T_ci(i),S(i),S_ci(i),m(i),m_ci(i));
    else
        line = sprintf('%d (%s) & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.2f\\pm%.2f$ \\\\\n',...
            i,seg_lbls{i==segments},u(i),u_ci(i),uc(i),uc_ci(i),T(i),T_ci(i),S(i),S_ci(i),m(i),m_ci(i));
    end
    line = strrep(line,'pm0.','pm.');
    line = strrep(line,'\','\\');
    fprintf(fid,['\t\t' line]);
end
fprintf(fid,'\t\t\\hline');

fclose(fid);





