% Script to generate segments for the glacier melt rate paper. This is a 
% summary of deployment times, locations, and depths for all deployments.
% The script generates everything within the 
% \begin{tabular} ... \end{tabular} commands.
%
% KJW
% 22 Oct 2026

clear

addpath('..')

% load data
dep_tbl = loadMSInfo(26:28);
seg_tbl = loadMSInfo(26:28,'segments');

% rows
rows = 1:22;
nrows = length(rows);

% column labels and units
col_lbls = {'Deployment','Segment','Location','Time','Duration','Depth'};
col_units = {'','','Lat,Lon (DD)','UTC','min','m'};

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
line1 = ['\t\t' line1(1:end-2) '\\\\\n'];
line2 = ['\t\t' line2(1:end-2) '\\\\\n'];

fid = fopen('../../../../../papers/glacier_2026/tables/table_segments_latex.txt','w');
fprintf(fid,line1);
fprintf(fid,line2);
fprintf(fid,'\t\t\\hline\n');

% time format
seg_tbl.Start.Format = 'yyyy-MM-dd HH:mm:SS';

% create rows
dep_prev = 0;
for i = 1:nrows
    dep_num = seg_tbl.Number(i);
    dep_row = dep_tbl.Number == dep_num;
    
    seg_num = seg_tbl.Window(i);
    lat = dep_tbl.Latitude(dep_row);
    lon = dep_tbl.Longitude(dep_row);
    ts = seg_tbl.Start(i);
    dur = round(minutes(seg_tbl.End(i)-ts),1);
    depth = round(seg_tbl.depth(i),1);

    if dep_num ~= dep_prev
        line = ['\t\t' sprintf('%d & %d & %.5f,%.5f & %s & %.1f & %.1f',dep_num-25,seg_num,lat,lon,ts,dur,depth)];
    else
        line = ['\t\t' sprintf(' & %d &  & %s & %.1f & %.1f',seg_num,ts,dur,depth)];
    end

    % print to file
    line = [line ' \\\\\n'];               
    fprintf(fid,line);

    dep_prev = dep_num;
end
fprintf(fid,'\t\t\\hline');

fclose(fid);





