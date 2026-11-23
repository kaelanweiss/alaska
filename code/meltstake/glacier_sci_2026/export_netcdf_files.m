% Script to output all of the processed meltstake data as netCDF files for
% Science 2026 submission.
%
% KJW
% 18 Nov 2026

clear

% load data
dep_tbl = loadMSInfo(26:28);

% data location
proc_dir = 'F:meltstake/data/proc';
raw_dir = 'F:meltstake/data/raw';
nc_dir = 'F:meltstake/data/netcdf/glacier_science_2026';

% load ctd
load(fullfile(proc_dir,'glacier_2024_ctd.mat'));

% bad salinity times
bad_sal_times = [datetime(2024,7,17,0,58,58.5) datetime(2024,7,17,1,1,7.5);...
                 datetime(2024,7,17,2,22,18.5) datetime(2024,7,17,2,22,36.5);...
                 datetime(2024,7,15,20,21,8) datetime(2024,7,15,20,35,10);...
                 datetime(2024,7,15,20,49,24) datetime(2024,7,15,20,56,40);...
                 ];

dep_lbls = [21 44 52];

% loop through deployments
ndeps = size(dep_tbl,1);
for i = 1:ndeps
    % load segment info
    seg_tbl = loadMSInfo(dep_tbl.Number(i),'segments');
    seg_tbl = seg_tbl(seg_tbl.include,:);

    % time limits
    t1 = seg_tbl.Start(1);
    t2 = seg_tbl.End(end);

    % adcp range info
    rmin = min(seg_tbl.rmax);
    rmean = mean(seg_tbl.rmax);

    % load data
    load(fullfile(raw_dir,dep_tbl.Folder{i},'adcp','adcp.mat'))
    load(fullfile(raw_dir,dep_tbl.Folder{i},'rbr','T.mat'))
    load(fullfile(raw_dir,dep_tbl.Folder{i},'hobo','hobo.mat'))
    load(fullfile(raw_dir,dep_tbl.Folder{i},'adv','adv.mat'))
    
    adcp.burst.processing.cor_min
    %%%% temperature %%%%
    time = T(1).time;
    idxt = time >= t1 & time <= t2;
    T_all = [T(1).values T(2).values T(3).values];
    % trim
    time = time(idxt);
    T_all = T_all(idxt,:);
    y_sensor = 0.01*[10 25 40]';
    % write file
    nc_name = fullfile(nc_dir,sprintf('XeitlSit_Jul2024_meltstake_deploy_%dm_temperature.nc',dep_lbls(i)));
    if exist(nc_name,'file')
        delete(nc_name)
    end
    % metadata
    writeVar(nc_name,'meltstake_depth',dep_lbls(i),{'meltstake_depth'},'m')
    writeVar(nc_name,'start_time',sprintf('%s UTC',t1),{'start_time'},'')
    writeVar(nc_name,'end_time',sprintf('%s UTC',t2),{'end_time'},'')
    % dimensions
    writeVar(nc_name,'time',posixtime(time),{'time'},'')
    ncwriteatt(nc_name,'time','format','posix')
    writeVar(nc_name,'y_sensor',y_sensor,{'y_sensor'},'m')
    % data
    writeVar(nc_name,'T',T_all,{'time','y_sensor'},'degC')

    %%%% hobo salinity %%%%
    hobo = hobo(2);
    time = hobo.time;
    % trim
    idx_bad_sal = false(size(time));
    for j = 1:size(bad_sal_times,1)
        tj1 = bad_sal_times(j,1);
        tj2 = bad_sal_times(j,2);
        idx_bad_sal = idx_bad_sal | (time>=tj1 & time<=tj2);
    end
    hobo.S(idx_bad_sal) = nan;

    idxt = time >= t1 & time <= t2;
    time = time(idxt);
    T_hobo = hobo.T_cal(idxt);
    C_hobo = hobo.C_cal(idxt);
    S_hobo = hobo.S(idxt);
    % write file
    nc_name = fullfile(nc_dir,sprintf('XeitlSit_Jul2024_meltstake_deploy_%dm_salinity.nc',dep_lbls(i)));
    if exist(nc_name,'file')
        delete(nc_name)
    end
    % metadata
    writeVar(nc_name,'meltstake_depth',dep_lbls(i),{'meltstake_depth'},'m')
    writeVar(nc_name,'start_time',sprintf('%s UTC',t1),{'start_time'},'')
    writeVar(nc_name,'end_time',sprintf('%s UTC',t2),{'end_time'},'')
    % dimensions
    writeVar(nc_name,'time',posixtime(time),{'time'},'')
    ncwriteatt(nc_name,'time','format','posix')
    % data
    writeVar(nc_name,'T',T_hobo,{'time'},'degC')
    writeVar(nc_name,'C',C_hobo,{'time'},'mS/cm')
    writeVar(nc_name,'S',S_hobo,{'time'},'psu')

    %%%% adcp %%%%
    time = adcp.burst.time;
    idxt = time >= t1 & time <= t2;
    idxr = adcp.burst.range <= rmin;
    % trim
    time = time(idxt);
    y = flip(rmean/cosd(25) - double(adcp.burst.range(idxr)))';

    uvw = flip(adcp.burst.vel_ice(idxt,idxr,[1 3 2]),2);
    pitch = interp1(adcp.attitude.time,adcp.attitude.pitch_corrected,time);
    roll = interp1(adcp.attitude.time,adcp.attitude.roll_corrected,time);
    pres = interp1(adcp.attitude.time,adcp.attitude.pressure,time);

    pitch(pitch<-180) = pitch(pitch<-180) + 360;
    roll(roll<-180) = roll(roll<-180) + 360;

    % write file
    nc_name = fullfile(nc_dir,sprintf('XeitlSit_Jul2024_meltstake_deploy_%dm_velocity.nc',dep_lbls(i)));
    if exist(nc_name,'file')
        delete(nc_name)
    end
    % metadata
    writeVar(nc_name,'meltstake_depth',dep_lbls(i),{'meltstake_depth'},'m')
    writeVar(nc_name,'start_time',sprintf('%s UTC',t1),{'start_time'},'')
    writeVar(nc_name,'end_time',sprintf('%s UTC',t2),{'end_time'},'')
    % dimensions
    writeVar(nc_name,'time',posixtime(time),{'time'},'')
    ncwriteatt(nc_name,'time','format','posix')
    writeVar(nc_name,'y',y,{'y'},'m')
    % data
    writeVar(nc_name,'u',uvw(:,:,1),{'time','y'},'m/s')
    writeVar(nc_name,'v',uvw(:,:,2),{'time','y'},'m/s')
    writeVar(nc_name,'w',uvw(:,:,3),{'time','y'},'m/s')
    writeVar(nc_name,'pressure',pres,{'time'},'dbar')
    writeVar(nc_name,'pitch',pitch,{'time'},'deg')
    % writeVar(nc_name,'roll',roll,{'time'},'deg')

    %%%% adv %%%%
    time = adv.pck_time;
    idxt = time >= t1 & time <= t2;
    % trim
    time = time(idxt);
    dist = adv.pck_dist(find(idxt,1),:)';
    amp = adv.pck_amp(idxt,:,:);

    % write file
    nc_name = fullfile(nc_dir,sprintf('XeitlSit_Jul2024_meltstake_deploy_%dm_backscatter.nc',dep_lbls(i)));
    if exist(nc_name,'file')
        delete(nc_name)
    end
    % metadata
    writeVar(nc_name,'meltstake_depth',dep_lbls(i),{'meltstake_depth'},'m')
    writeVar(nc_name,'start_time',sprintf('%s UTC',t1),{'start_time'},'')
    writeVar(nc_name,'end_time',sprintf('%s UTC',t2),{'end_time'},'')
    % dimensions
    writeVar(nc_name,'time',posixtime(time),{'time'},'')
    ncwriteatt(nc_name,'time','format','posix')
    writeVar(nc_name,'dist',dist,{'dist'},'mm')
    % data
    writeVar(nc_name,'amp1',amp(:,:,1),{'time','dist'},'count')
    writeVar(nc_name,'amp2',amp(:,:,2),{'time','dist'},'count')
    writeVar(nc_name,'amp3',amp(:,:,3),{'time','dist'},'count')

    %%%% CTD %%%%
    nz = length(ctd(i).depth_mean);
    T = ctd(i).T(1:nz,:);
    S = ctd(i).S(1:nz,:);
    rho = ctd(i).PD(1:nz,:)+1000;
    N = ctd(i).N;

    % write file
    nc_name = fullfile(nc_dir,sprintf('XeitlSit_Jul2024_meltstake_deploy_%dm_fjordCTD.nc',dep_lbls(i)));
    if exist(nc_name,'file')
        delete(nc_name)
    end
    % metadata
    writeVar(nc_name,'meltstake_depth',dep_lbls(i),{'meltstake_depth'},'m')
    writeVar(nc_name,'start_time',sprintf('%s UTC',t1),{'start_time'},'')
    writeVar(nc_name,'end_time',sprintf('%s UTC',t2),{'end_time'},'')
    % dimensions
    writeVar(nc_name,'time',posixtime(ctd(i).time'),{'time'},'')
    ncwriteatt(nc_name,'time','format','posix')
    writeVar(nc_name,'depth',ctd(i).depth_mean,{'depth'},'m')
    % data
    writeVar(nc_name,'lat',ctd(i).lat',{'time'},'dec deg')
    writeVar(nc_name,'lon',ctd(i).lon',{'time'},'dec deg')
    writeVar(nc_name,'T',T,{'depth','time'},'degC')
    writeVar(nc_name,'S',S,{'depth','time'},'psu')
    writeVar(nc_name,'rho_pot',rho,{'depth','time'},'kg/m^3')
    writeVar(nc_name,'N',N,{'depth'},'cps')

end







% %%%%%%%%
% 
% % collect folder names and number of segments in each deployment
% folders = unique(ms_tbl.Folder,'stable');
% n_segs = nan(length(folders),1);
% for i = 1:length(folders)
%     n_segs(i) = sum(strcmp(ms_tbl.Folder,folders{i}));
% end
% 
% % loop through deployment folders
% seg = 1; % keep track of segment number
% for i = 1:length(folders)
%     dep_dir = fullfile(proc_dir,folders{i});
% 
%     % load adv data
%     load(fullfile(raw_dir,folders{i},'adv','adv.mat'),'adv')
% 
%     % loop through segments
%     for j = 1:n_segs(i)
%         % grab the correct pck data
%         t1 = ms_tbl.Start(seg);
%         t2 = ms_tbl.End(seg);
%         idxt = adv.pck_time >= t1 & adv.pck_time <= t2;
%         time_echo = adv.pck_time(idxt);
%         dist_echo = adv.pck_dist(find(idxt,1),:)';
%         amp_echo = adv.pck_amp(idxt,:,:);
% 
%         % load the stuff
%         load(fullfile(dep_dir,sprintf('ctd%d.mat',j)))
%         load(fullfile(dep_dir,sprintf('T%d.mat',j)))
%         load(fullfile(dep_dir,sprintf('adcp%d.mat',j)))
% 
%         % segment depth
%         seg_depth = ms_tbl.depth(seg);
% 
%         % define dimensions
%         time_T = T.time;
%         sensor_number_T = (1:size(T.T,2))';
%         time_S = ctd.t0';
%         z_S = -ctd.z;
%         time_uvw = adcp.burst.time;
%         y_uvw = ms_tbl.rmax(seg)/cosd(25) - adcp.burst.range';
% 
%         % change date format
%         t_fmt = 'uuuu-MM-dd HH:mm:ss';
%         t1.Format = t_fmt;
%         t2.Format = t_fmt;
% 
%         % convert datetimes to epochtime
%         time_T = posixtime(time_T);
%         time_S = posixtime(time_S);
%         time_uvw = posixtime(time_uvw);
%         time_echo = posixtime(time_echo);
% 
%         % write file
%         nc_name = fullfile(nc_dir,sprintf('XeitlSit_2023_meltstake_segment_%02d.nc',seg));
%         if exist(nc_name,'file')
%             delete(nc_name)
%         end
%         % metadata
%         writeVar(nc_name,'segment_number',seg,{'segment_number'},'')
%         writeVar(nc_name,'meltstake_depth',seg_depth,{'meltstake_depth'},'m')
%         writeVar(nc_name,'start_time',sprintf('%s UTC',t1),{'start_time'},'')
%         writeVar(nc_name,'end_time',sprintf('%s UTC',t2),{'end_time'},'')
%         % dimensions
%         writeVar(nc_name,'time_T',time_T,{'time_T'},'')
%         writeVar(nc_name,'sensor_number_T',sensor_number_T,{'sensor_number_T'},'')
%         writeVar(nc_name,'time_S',time_S,{'time_S'},'')
%         writeVar(nc_name,'z_S',z_S,{'z_S'},'m')
%         writeVar(nc_name,'time_uvw',time_uvw,{'time_uvw'},'')
%         writeVar(nc_name,'y_uvw',y_uvw,{'y_uvw'},'m')
%         writeVar(nc_name,'time_echo',time_echo,{'time_echo'},'')
%         writeVar(nc_name,'dist_echo',dist_echo,{'dist_echo'},'mm')
%         % T
%         writeVar(nc_name,'T',T.T,{'time_T','sensor_number_T'},'degC')
%         % S
%         writeVar(nc_name,'S',ctd.S',{'time_S','z_S'},'psu')
%         % vel
%         writeVar(nc_name,'u',adcp.burst.vel_ice(:,:,1),{'time_uvw','y_uvw'},'m/s')
%         writeVar(nc_name,'v',adcp.burst.vel_ice(:,:,3),{'time_uvw','y_uvw'},'m/s')
%         writeVar(nc_name,'w',adcp.burst.vel_ice(:,:,2),{'time_uvw','y_uvw'},'m/s')
%         % echo
%         writeVar(nc_name,'echo1',amp_echo(:,:,1),{'time_echo','dist_echo'},'counts')
%         writeVar(nc_name,'echo2',amp_echo(:,:,2),{'time_echo','dist_echo'},'counts')
%         writeVar(nc_name,'echo3',amp_echo(:,:,3),{'time_echo','dist_echo'},'counts')
% 
%         % add format to time dimensions (unnecessary with posixtime)
%         time_vars = {'time_T','time_S','time_uvw','time_echo'};
%         for k = 1:length(time_vars)
%             ncwriteatt(nc_name,time_vars{k},'format','posix')
%         end
% 
%         % next segment
%         seg = seg + 1;
%     end
% end

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