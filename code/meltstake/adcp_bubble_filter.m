% Script to clean up bubbles/sediment from ADCP records.
%
% KJW
% 27 Mar 2025

clear

proc_dir = 'F:/meltstake/data/proc';

ms_tbl = loadMSInfo(26:28,'manualwindows');

[dep_nums,uidx] = unique(ms_tbl.Number);
dep_names = ms_tbl.Folder(uidx);
ndeps = length(dep_nums);

% filter parameters
nsigma = 3;
npass = 2;
width_T = 5; %s
ovrlp_T = width_T/3; %s
detrend_flag = true;

% loop through deployments
for i = 1%:ndeps
    seg_nums = ms_tbl.Window(ms_tbl.Number==dep_nums(i));
    
    % loop through segments
    for j = 1:length(seg_nums)
        % load
        load(fullfile(proc_dir,dep_names{i},sprintf('adcp%d.mat',j)))
        % plot data
        ax = adcpQuickPlot(figure(1),adcp,'vel_unw',max(abs(extrema(adcp.burst.vel)))*[-1 1],NaT,[0 1],1);
        
        % preallocate bad indices
        sidx1 = false(size(adcp.burst.vel));
        sidx2 = false(size(adcp.burst.vel));
        
        % sigma filter along time axis
        % calculate filter width and overlap based on sample rate
        width = round(width_T*adcp.burst.samplerate);
        ovrlp = round(ovrlp_T*adcp.burst.samplerate);
        % loop through beams
        for k = 1:size(adcp.burst.vel,3)
            [~,sidxk] = runningSigmaFilter(adcp.burst.vel_unw(:,:,k),nsigma,npass,width,ovrlp,detrend_flag);
            sidx1(:,:,k) = sidxk;
            sum(sidxk,'all')/numel(sidxk)
        end

        % sigma filter over profiles
        % beams
        for k = 1:size(adcp.burst.vel,3)
            % time
            for m = 1:size(adcp.burst.vel,1)
                vel_prof = adcp.burst.vel_unw(m,:,k);
                vel_mean = hannFilter(vel_prof,3);
                vel_peaks = vel_prof - vel_mean;
                [~,sidxmk] = sigmaFilter(vel_peaks,2,2,1);
                sidx2(m,:,k) = sidxmk;
            end
        end

        % combine
        sidx = sidx1 & sidx2;
        
        % add bad cells to plots
        for k = 1:size(adcp.burst.vel,3)
            hold(ax(k),'on')
            % loop through cells
            for m = 1:size(adcp.burst.vel,2)
                ym = adcp.burst.range(m);
                tm = adcp.burst.time(sidx(:,m,k));
                plot(ax(k),tm+seconds(0.5/adcp.burst.samplerate),ym*ones(size(tm))+adcp.burst.cellsize/2,'o','color',colors(4))
            end
        end
        
        a = 0;
        
    end
end