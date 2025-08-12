function tlim = parseDeploymentTimes(depname)
% Function to parse time limits of a deployment in a "deployment_times.txt"
% file.
%
% tlim = parseDeploymentTimes(depname)
%
% Input
%   depname: folder containing "deployment_times.txt" file
%
% Output
%   tlim: 1x2 datenum vector of beginning and end of deployment
%
% KJW
% 16 Sep 2022

tlim = [nan nan];
fid = fopen(fullfile(depname,'deployment_times.txt'),'r');
for i = 1:2
    line = fgetl(fid);
    while line(1) == '#'
        line = fgetl(fid);
    end
    vals = strsplit(line,',');
    tlim(i) = datenum(strjoin(vals(1:6),'-'),'yyyy-mm-dd-HH-MM-SS');
end
fclose(fid);