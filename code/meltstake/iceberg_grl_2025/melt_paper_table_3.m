% Script to generate table 2 for the melt rate paper. This is a summary of
% deployment times and depth for all deployments. The script generates
% everything within the \begin{tabular} ... \end{tabular} commands.
%
% KJW
% 10 Apr 2024

clear

addpath('..')

% load data
ms_tbl = loadMSInfo(1:19,'manualwindows');

% rows
rows = 1:31;
nrows = length(rows);
n_splt = ceil(nrows/2);


% column labels and units
col_lbls = {'Segment','Time','Duration','Depth'};
col_units = {'','UTC','min','m'};

% header lines
ncols = length(col_lbls);
line1 = '';
line2 = '';
for i = 1:ncols
    line1 = [line1 sprintf('%s & ',col_lbls{i})];
    if ~isempty(col_units{i})
        line2 = [line2 sprintf('[%s] & ',col_units{i})];
    else
        line2 = [line2 ' & '];
    end
end
line1 = ['\t\t' line1 line1(1:end-2) '\\\\\n'];
line2 = ['\t\t' line2 line2(1:end-2) '\\\\\n'];

fid = fopen('../../../../../papers/weiss_et_al_2024/table_deployments_latex.txt','w');
fprintf(fid,line1);
fprintf(fid,line2);
fprintf(fid,'\t\t\\hline\n');

% time format
ms_tbl.Start.Format = 'yyyy/MM/dd HH:mm';

% cases
segments = [8 3 29 9];
seg_lbls = upper({'a','b','c','d'});

% create rows
for i = 1:n_splt
    line = '\t\t';
    for j = 1:2
        idx = n_splt*(j-1)+i;
        tstamp = ms_tbl.Start(idx);
        if idx <= nrows
            if ~any(idx==segments)
                line = [line sprintf('%d & %s & %d & %.1f',idx,tstamp,...
                    round(ms_tbl.Duration(idx)),ms_tbl.depth(idx))];
            else
                line = [line sprintf('%d (%s) & %s & %d & %.1f',idx,seg_lbls{idx==segments},tstamp,...
                    round(ms_tbl.Duration(idx)),ms_tbl.depth(idx))];
            end
        end
        if j == 1
            line = [line ' & '];
        end
    end
    
    line = [line ' \\\\\n'];               
%     line = strrep(line,'\','\\');
    fprintf(fid,line);
end
fprintf(fid,'\t\t\\hline');

fclose(fid);





