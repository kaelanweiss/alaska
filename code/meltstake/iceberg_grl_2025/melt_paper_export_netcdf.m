% Script to output all of the processed meltstake data as netCDF files for
% GRL submission.
%
% KJW
% 16 Jan 2025

clear

% load data
ms_tbl = loadMSInfo(1:19,'manualwindows');

% data location
proc_dir = 'F:meltstake/data/proc';
raw_dir = 'F:meltstake/data/raw';
nc_dir = 'F:meltstake/data/netcdf';

% collect folder names and number of segments in each deployment
folders = unique(ms_tbl.Folder,'stable');
n_segs = nan(length(folders),1);
for i = 1:length(folders)
    n_segs(i) = sum(strcmp(ms_tbl.Folder,folders{i}));
end

% loop through deployment folders
seg = 1; % keep track of segment number
for i = 1:length(folders)
    dep_dir = fullfile(proc_dir,folders{i});

    % load adv data
    load(fullfile(raw_dir,folders{i},'adv','adv.mat'),'adv')

    % loop through segments
    for j = 1:n_segs(i)
        % grab the correct pck data
        t1 = ms_tbl.Start(seg);
        t2 = ms_tbl.End(seg);
        idxt = adv.pck_time >= t1 & adv.pck_time <= t2;
        time_echo = adv.pck_time(idxt);
        dist_echo = adv.pck_dist(find(idxt,1),:)';
        amp_echo = adv.pck_amp(idxt,:,:);

        % load the stuff
        load(fullfile(dep_dir,sprintf('ctd%d.mat',j)))
        load(fullfile(dep_dir,sprintf('T%d.mat',j)))
        load(fullfile(dep_dir,sprintf('adcp%d.mat',j)))

        % segment depth
        seg_depth = ms_tbl.depth(seg);

        % define dimensions
        time_T = T.time;
        sensor_number_T = (1:size(T.T,2))';
        time_S = ctd.t0';
        z_S = -ctd.z;
        time_uvw = adcp.burst.time;
        y_uvw = ms_tbl.rmax(seg)/cosd(25) - adcp.burst.range';

        % change date format
        t_fmt = 'uuuu-MM-dd HH:mm:ss';
        t1.Format = t_fmt;
        t2.Format = t_fmt;

        % convert datetimes to epochtime
        time_T = posixtime(time_T);
        time_S = posixtime(time_S);
        time_uvw = posixtime(time_uvw);
        time_echo = posixtime(time_echo);

        % write file
        nc_name = fullfile(nc_dir,sprintf('XeitlSit_2023_meltstake_segment_%02d.nc',seg));
        if exist(nc_name,'file')
            delete(nc_name)
        end
        % metadata
        writeVar(nc_name,'segment_number',seg,{'segment_number'},'')
        writeVar(nc_name,'meltstake_depth',seg_depth,{'meltstake_depth'},'m')
        writeVar(nc_name,'start_time',sprintf('%s UTC',t1),{'start_time'},'')
        writeVar(nc_name,'end_time',sprintf('%s UTC',t2),{'end_time'},'')
        % dimensions
        writeVar(nc_name,'time_T',time_T,{'time_T'},'')
        writeVar(nc_name,'sensor_number_T',sensor_number_T,{'sensor_number_T'},'')
        writeVar(nc_name,'time_S',time_S,{'time_S'},'')
        writeVar(nc_name,'z_S',z_S,{'z_S'},'m')
        writeVar(nc_name,'time_uvw',time_uvw,{'time_uvw'},'')
        writeVar(nc_name,'y_uvw',y_uvw,{'y_uvw'},'m')
        writeVar(nc_name,'time_echo',time_echo,{'time_echo'},'')
        writeVar(nc_name,'dist_echo',dist_echo,{'dist_echo'},'mm')
        % T
        writeVar(nc_name,'T',T.T,{'time_T','sensor_number_T'},'degC')
        % S
        writeVar(nc_name,'S',ctd.S',{'time_S','z_S'},'psu')
        % vel
        writeVar(nc_name,'u',adcp.burst.vel_ice(:,:,1),{'time_uvw','y_uvw'},'m/s')
        writeVar(nc_name,'v',adcp.burst.vel_ice(:,:,3),{'time_uvw','y_uvw'},'m/s')
        writeVar(nc_name,'w',adcp.burst.vel_ice(:,:,2),{'time_uvw','y_uvw'},'m/s')
        % echo
        writeVar(nc_name,'echo1',amp_echo(:,:,1),{'time_echo','dist_echo'},'counts')
        writeVar(nc_name,'echo2',amp_echo(:,:,2),{'time_echo','dist_echo'},'counts')
        writeVar(nc_name,'echo3',amp_echo(:,:,3),{'time_echo','dist_echo'},'counts')

        % add format to time dimensions (unnecessary with posixtime)
        time_vars = {'time_T','time_S','time_uvw','time_echo'};
        for k = 1:length(time_vars)
            ncwriteatt(nc_name,time_vars{k},'format','posix')
        end

        % next segment
        seg = seg + 1;
    end
end

%%%%% subfunctions %%%%%
function writeVar(fname,name,data,dim_names,units)
    % build dimensions cell array
    dims = cell(1,2*length(dim_names));
    for i = 1:length(dim_names)
        dims{2*i-1} = dim_names{i};
        dims{2*i} = size(data,i);
    end
    % create file/variable
    if isstr(data) || ischar(data)
        nccreate(fname,name,'dimensions',dims,'datatype','string','format','netcdf4')
    else
        nccreate(fname,name,'dimensions',dims,'format','netcdf4')
    end
    % write data
    ncwrite(fname,name,data)
    % write units
    if ~strcmp(units,'')
        ncwriteatt(fname,name,'units',units)
    end

end