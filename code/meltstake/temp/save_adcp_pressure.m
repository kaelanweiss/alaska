% Script to extract and save ADCP pressure data in a more generalized
% format
%
% KJW
% 1 Jul 2025

clear

raw_dir = 'F:/meltstake/data/raw';

ms_tbl = loadMSInfo(26:28);
dep_names = ms_tbl.Folder;

for i = 1:length(dep_names)
    % load adcp data
    load(fullfile(raw_dir,dep_names{i},'adcp','adcp.mat'))
    time = adcp.attitude.time;
    pres = adcp.attitude.pressure;
    roll = -adcp.attitude.pitch; % when +y_inst points down, "pitch" is zero and increases with negative rotation about the +z_inst axis
    pitch = (adcp.attitude.roll-270); % "roll" is positive rotation about +x_inst, at 270 when +y_inst points down
    pitch(pitch<-180) = pitch(pitch<-180)+360;
    n = length(time);
    
    % reformat time string
    time.Format = 'uuuu-MM-dd''T''HH:mm:ss.SSS';
    
    % create csv file
    fid = fopen(fullfile('temp',[dep_names{i} '.csv']),'w');
    fprintf(fid,'timestamp,pressure,pitch,roll\n');
    fprintf(fid,'UTC,dbar,deg,deg\n');

    fprintf('Writing %s.csv...',dep_names{i})

    for j = 1:n
        fprintf(fid,'%s,%.3f,%.1f,%.1f\n',time(j),round(pres(j),3),round(pitch(j),1),round(roll(j),1));
    end

    fclose(fid);

    fprintf('done.\n')

end