% Script to generate results table for the glacier melt rate paper. This is
% a summary of deployment data for all deployments, including predicted
% melt rate. The script generates everything within the
% \begin{tabular} ... \end{tabular} commands.
%
% KJW
% 22 Oct 2026

clear

addpath('..')

% load data
ms_tbl = loadMSInfo(26:28,'segments');

% extract fields
m = msTable2Vector(ms_tbl.m);
[m_ci_l,m_ci_h] = msTable2Vector(ms_tbl.m_ci);
[T,T_ci] = msTable2Vector(ms_tbl.T);
[S,S_ci] = msTable2Vector(ms_tbl.S);
[u,u_ci] = msTable2Vector(ms_tbl.u_max);
uc = msTable2Vector(ms_tbl.u_coarse);
m_ci = [m_ci_l,m_ci_h];

% load melt rates
mod_glac = load('../../../data/melt_models_glacier_hann_1sec.mat');
m3eqn = mod_glac.m_mean(:,10);
m3eqn_coarse = mod_glac.m_mean(:,11);
Tb = mod_glac.Tb_mean;
Sb = mod_glac.Sb_mean;

% calculate Tf
Tf = T-Tb;

% convert melt rate to m/day
m = m*.24;
m_ci = m_ci*.24;

% calculate melt rate [low high]
m_low = m-m_ci(:,1);
m_high = m+m_ci(:,2);

% round to desired number of decimal places
m = round(m,2);
m_low = round(m_low,2);
m_high = round(m_high,2);
Tf = round(Tf,1);
S = round(S,1);
u = round(u,3);
uc = round(uc,3);
Tb = round(Tb,2);
Sb = round(Sb,1);
m3eqn = round(m3eqn,2);
m3eqn_coarse = round(m3eqn_coarse,2);

% rows
rows = 1:22;
nrows = length(rows);

% column labels and units
col_lbls = {'Seg','$\\overline{U_\\mathrm{inst}}$','$\\overline{U_\\mathrm{mean}}$','$\\overline{T-T_b}$','$\\overline{S}$',...
            '$m_\\mathrm{obs}$','$\\overline{T_b}$','$\\overline{S_b}$','$m_\\mathrm{3eqn}$','$m_\\mathrm{3eqn,coarse}$'};
col_units = {'','m/s','m/s','K','psu','m/day','$^\\circ$C','psu','m/day','m/day'};

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

fid = fopen('../../../../../papers/glacier_2026/tables/table_results_latex.txt','w');
fprintf(fid,line1);
fprintf(fid,line2);
fprintf(fid,'\t\t\\hline\n');

% create rows
for i = 1:nrows
    dep_num = ms_tbl.Number(i)-25;
    seg_num = ms_tbl.Window(i);

    % if seg_num > 1
    %     dep_num = ' ';
    % end
    
    line = sprintf('%d.%d & %.3f & %.3f & %.1f & %.1f & %.2f (%.2f, %.2f) & %.2f & %.1f & %.2f & %.2f',...
                   dep_num,seg_num,u(i),uc(i),Tf(i),S(i),m(i),m_low(i),m_high(i),Tb(i),Sb(i),m3eqn(i),m3eqn_coarse(i));

    fprintf(fid,['\t\t' line ' \\\\\n']);
end

% % create rows
% for i = 1:nrows
%     if ~any(i==segments)
%         line = sprintf('%d & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.2f\\pm%.2f$ \\\\\n',...
%             i,u(i),u_ci(i),uc(i),uc_ci(i),T(i),T_ci(i),S(i),S_ci(i),m(i),m_ci(i));
%     else
%         line = sprintf('%d (%s) & $%.3f\\pm%.3f$ & $%.3f\\pm%.3f$ & $%.1f\\pm%.1f$ & $%.1f\\pm%.1f$ & $%.2f\\pm%.2f$ \\\\\n',...
%             i,seg_lbls{i==segments},u(i),u_ci(i),uc(i),uc_ci(i),T(i),T_ci(i),S(i),S_ci(i),m(i),m_ci(i));
%     end
%     line = strrep(line,'pm0.','pm.');
%     line = strrep(line,'\','\\');
%     fprintf(fid,['\t\t' line]);
% end
fprintf(fid,'\t\t\\hline');

fclose(fid);





