% Script to generate melt rate correction table for the glacier melt rate 
% paper. This is a summary of the inputs and resulting melt bias for all 
% deployments. The script generates everything within the
% \begin{tabular} ... \end{tabular} commands.
%
% KJW
% 18 Nov 2026

clear

addpath('..')

% load data
ms_tbl = loadMSInfo(26:28,'pitch_correction');

% extract fields
p = msTable2Vector(ms_tbl.pitch_mean);
pr = msTable2Vector(ms_tbl.pitch_rate);
beta = msTable2Vector(ms_tbl.ice_slope);
mb = msTable2Vector(ms_tbl.melt_bias);

% convert melt rate to m/day
mb = mb*.24;

% round to desired number of decimal places
p = round(p,2);
pr = round(pr,2);
beta = round(beta,1);
mb = round(mb,2);

% rows
rows = 1:22;
nrows = length(rows);

% column labels and units
col_lbls = {'Segment','$\\phi$','$\\dot{\\phi}$','$\\beta$','$\\Delta m$'};
col_units = {'','deg','deg\\per{hr}','deg','m\\per{day}'};

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

fid = fopen('../../../../../papers/glacier_2026/tables/table_melt_correction_latex.txt','w');
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
    
    line = sprintf('%d.%d & %.2f & %.2f & %.1f & %.2f',...
                   dep_num,seg_num,p(i),pr(i),beta(i),mb(i));

    fprintf(fid,['\t\t' line ' \\\\\n']);
end

fprintf(fid,'\t\t\\hline');

fclose(fid);





