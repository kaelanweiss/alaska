% script to calculate an outer velocity scale from meltstake ADCP
% observations to use as input to the melt parameterization
%
% 1) take a swath of velocity data a certain distance from the ADCP
% 2) average bins together to get a single time series of u, v, w
% 3) add in quadrature (u,v,w) and (u,w)
% 4) make some plots
% ...
%
% KJW
% 15 Aug 2023

clear

%%% MAy 2023 %%%
% proc_path = 'F:/AK_202305/adcp/proc';
% 
% ndeps = [3];
% 
% max_ranges = {[0.45 0.5 0.55]};
% 
% time_ranges = {{{{'0105','0200'}},{{'2040','2136'},{'2141','2159'}},{{'0147','0220'}}}};


%%% July 2023 %%%
% proc_path = 'F:/AK_202307/adcp/proc';
% 
% ndeps = [1 4 5];
% 
% max_ranges = {[0.65],...
%               [0.6, 0.4, 0.4, 0.3],...
%               [0.6, 0.3, 0.55, 0.55, 0.55]};
% 
% time_ranges = {{{{''}}},...
%                {{{'1233','1415'}},{{'1630','1645'}},{{'1130','1330'}},{{''}}},...
%                {{{'1445','1510'}},{{'1712','1752'}},{{'1400','1520'},{'1530','1700'}},{{'1100','1230'}},{{'1530','1645'},{'1645','1740'}}}};

%%% September 2023 %%%
proc_path = 'F:/AK_202309/data/meltstake/proc';

ndeps = [6];

max_ranges = {[]};

time_ranges = {{{{'0105','0200'}},{{'2040','2136'},{'2141','2159'}},{{'0147','0220'}}}};



u_scale = time_ranges;
RLIM = 0.3;

for i = 1:length(ndeps) % loop through meltstakes
    for j = 1:ndeps(i) % loop through deployments
        % load data
        %load(fullfile(proc_path,sprintf('ms%02d_%02d_adcp.mat',i,j)))
        load(fullfile(proc_path,sprintf('dep%d_adcp.mat',j)))
        
        % calculate mean
        idxr = (adcp.burst.range-1e-3)<=RLIM; % pad to be inclusive
        vel = squeeze(mean(adcp.burst.vel_ice(:,idxr,1:3),2,'omitnan'));

        % calculate magnitude
        vel_mag = nan(size(vel,1),2);
        vel_mag(:,1) = sqrt(sum(vel.^2,2));
        vel_mag(:,2) = sqrt(sum(vel(:,[1 2]).^2,2));

        % calculate velocity scales
        day = datestr(adcp.burst.time(1),'yyyymmdd');
        for k = 1:length(time_ranges{i}{j})
            % find slice
            trk = time_ranges{i}{j}{k};
            t1 = datenum(0);
            t2 = datenum(Inf);
            if ~isempty(trk{1})
                t1 = datenum([day trk{1}],'yyyymmddHHMM');
                t2 = datenum([day trk{2}],'yyyymmddHHMM');
            end
            idxk = adcp.burst.time>=t1 & adcp.burst.time<=t2;

            % calculate scale
            u_scale{i}{j}{k} = [nan nan nan nan];
            u_scale{i}{j}{k}(1) = mean(vel_mag(idxk,1),'omitnan');
            u_scale{i}{j}{k}(2) = mean(vel_mag(idxk,2),'omitnan');
            u_scale{i}{j}{k}(3) = sqrt(mean(vel_mag(idxk,1).^2,'omitnan'));
            u_scale{i}{j}{k}(4) = sqrt(mean(vel_mag(idxk,2).^2,'omitnan'));
        end



        % plot
        t = adcp.burst.time;
        figure(j); clf; hold on
        plot(t,vel)
        plot(t,vel_mag(:,1),'k-')
        plot(t,vel_mag(:,2),'k--')
        xlim(extrema(t))
        datetick('x','mmmdd HH:MM','keeplimits')
        legend(cat(2,adcp.burst.processing.ice_lbls(1:3),{'sqrt(u^2+v^2+w^2)','sqrt(u^2+w^2)'}))
        grid on
        box on
        drawnow
        pause(0.01)

    end
end

for i = 1:length(ndeps)
    fprintf('ms%02d\n',i)
    for j = 1:ndeps(i)
        fprintf('\tdep%02d\n',j)
        for k = 1:length(u_scale{i}{j})
            fprintf('\t\tsection%02d: %.3f m/s\n',k,u_scale{i}{j}{k}(2))
        end
    end
end


