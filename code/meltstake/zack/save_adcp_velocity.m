% Script to extract and save y-averaged ADCP velocity data in a more
% generalized format.
%
% KJW
% 23 Sep 2025

clear
addpath('..')

raw_dir = 'F:/meltstake/data/raw';

ms_tbl = loadMSInfo(26:28);
dep_names = ms_tbl.Folder;

rmax = 0.4;

for i = 1:length(dep_names)
    % load adcp data
    load(fullfile(raw_dir,dep_names{i},'adcp','adcp.mat'))
    time = adcp.burst.time;
    idxr = adcp.burst.range<=rmax;
    
    vel_ice = squeeze(mean(adcp.burst.vel_ice(:,idxr,[1 3 2]),2,'omitnan'));
    n = length(time);
    
    % reformat time string
    time.Format = 'uuuu-MM-dd''T''HH:mm:ss.SSS';
    
    % create csv file
    fid = fopen(fullfile('meltstake_mean_velocity',[dep_names{i} '_uvw.csv']),'w');
    fprintf(fid,'timestamp,u,v,w\n');
    fprintf(fid,'UTC,m/s,m/s,m/s\n');

    fprintf('Writing %s_uvw.csv...',dep_names{i})

    for j = 1:n
        fprintf(fid,'%s,%.4f,%.4f,%.4f\n',time(j),round(vel_ice(j,1),4),round(vel_ice(j,2),4),round(vel_ice(j,3),4));
    end

    fclose(fid);

    fprintf('done.\n')

end