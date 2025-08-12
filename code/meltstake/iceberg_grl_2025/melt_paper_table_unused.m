% Script to generate table 1 for the melt rate paper. This is a summary of
% a few meltstake sections. The script generates everything within the
% \begin{tabular} ... \end{tabular} commands.
%
% KJW
% 10 Apr 2024

clear

% load data
tbl_path = 'G:Shared drives\Ice-ocean-interactions\science\Grad Students\Kaelan\meltstake_deployments.xlsx';
ms_tbl = readtable(tbl_path,'sheet','manualwindows');
load ../../data/melt_models.mat

% cm/hr to m/dy
cmh2md = 0.24;

% sections
sections = [8 3 29 9];
section_lbls = {'crossflow','plume','wave','intermittent'};
nrows = length(sections);

% data
m_obs = msTable2Vector(ms_tbl.m(sections));
m_obs = m_obs*cmh2md;
m_obs_ci = ms_tbl.m_ci(sections);
m_obs_ci = m_obs_ci*cmh2md;

m_3eqn = m_mean(sections,1);
m_err = m_obs./m_3eqn;
m_err_ci = m_obs_ci./m_3eqn;

[T,T_ci] = msTable2Vector(ms_tbl.T(sections));
[u,u_ci] = msTable2Vector(ms_tbl.u0(sections));

% column labels and units
col_lbls = {'Case','$T$','$U$','$m_{obs}$','$m_{3eqn}$','$r_{3eqn}$'};
col_units = {'','$^\\circ$C','m/s','m/day','m/day',''};

ncols = length(col_lbls);

line1 = '\t\t';
line2 = '\t\t';
for i = 1:ncols
    line1 = [line1 sprintf('%s & ',col_lbls{i})];
    if ~isempty(col_units{i})
        line2 = [line2 sprintf('(%s) & ',col_units{i})];
    else
        line2 = [line2 ' & '];
    end
end
line1 = [line1(1:end-2) '\\\\\n'];
line2 = [line2(1:end-2) '\\\\\n'];

fid = fopen('../../../../papers/weiss_et_al_2024/table1_latex.txt','w');
fprintf(fid,line1);
fprintf(fid,line2);

% create rows
fprintf(fid,'\t\t\\hline\n');
for i = 1:nrows
    line = sprintf('\t\t%c (%s) & $%.1f\\pm%.1f$ & $%.3f\\pm%.3f$ & $%.2f\\pm%.2f$ & $%.2f$ & $%.1f\\pm%.1f$ \\\\\n',...
                   char(64+i),section_lbls{i},T(i),T_ci(i),u(i),u_ci(i),m_obs(i),m_obs_ci(i),m_3eqn(i),m_err(i),m_err_ci(i));
    line = strrep(line,'pm0.','pm.');
    line = strrep(line,'\','\\');
    fprintf(fid,line);
end
fprintf(fid,'\t\t\\hline');

fclose(fid);





