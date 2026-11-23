% Script to check the data exported to .nc for the glacier Science paper.
%
% KJW
% 19 Nov 2026

clear

% load info
nc_dir = 'F:meltstake/data/netcdf/glacier_science_2026';
nc_pre = 'XeitlSit_Jul2024_meltstake_deploy_';
flds = {'backscatter','fjordCTD','salinity','temperature','velocity'};
deps = {'21','44','52'};

%% load data into one big structure vector
clear data
data(3) = struct();
for i = 1:3
    fprintf('--- Deploy %s ---\n',deps{i})
    for fld = flds
        fprintf('%s\n',fld{1})
        data(i).(fld{1}) = struct();
        filej = fullfile(nc_dir,sprintf('%s%sm_%s.nc',nc_pre,deps{i},fld{1}));
        finfo = ncinfo(filej);
        var_names = {finfo.Variables.Name};
        for vark = var_names
            data(i).(fld{1}).(vark{1}) = ncread(filej,vark{1});
            if strcmp(vark{1},'time')
                data(i).(fld{1}).(vark{1}) = datetime(data(i).(fld{1}).(vark{1}),'convertfrom','posix');
            end
            if isempty(finfo.Variables(strcmp(vark,var_names)).Attributes)
                units = '';
            else
                units = finfo.Variables(strcmp(vark,var_names)).Attributes.Value;
            end
            fprintf('  %s (%s)\n',vark{1},units)
        end
    end
end

%% plot
for i = 1:3
    % backscatter
    figure(10*i+1); clf
    for j = 1:3
        subplot(1,3,j)
        pcolor(data(i).backscatter.time,data(i).backscatter.dist,data(i).backscatter.(['amp' num2str(j)])')
        shading flat
    end

    % fjordCTD
    figure(10*i+2); clf
    subplot(1,4,1)
    plot(data(i).fjordCTD.T,data(i).fjordCTD.depth)
    axis ij
    subplot(1,4,2)
    plot(data(i).fjordCTD.S,data(i).fjordCTD.depth)
    axis ij
    subplot(1,4,3)
    plot(data(i).fjordCTD.rho_pot,data(i).fjordCTD.depth)
    axis ij
    subplot(1,4,4)
    plot(data(i).fjordCTD.N,data(i).fjordCTD.depth)
    axis ij

    % salinity
    figure(10*i+3); clf
    subplot(3,1,1)
    plot(data(i).salinity.time,data(i).salinity.C)
    subplot(3,1,2)
    plot(data(i).salinity.time,data(i).salinity.T)
    subplot(3,1,3)
    plot(data(i).salinity.time,data(i).salinity.S)

    % temperature
    figure(10*i+4); clf
    hold on
    h_lbls = {};
    for j = 1:3
        h(j) = plot(data(i).temperature.time,data(i).temperature.T(:,j));
        h_lbls{j} = [num2str(round(100*data(i).temperature.y_sensor(j))) ' cm'];
    end
    legend(h,h_lbls)

    % velocity
    figure(10*i+5); clf
    adcp = struct('burst',struct());
    adcp.burst.samplerate = round(1/seconds(diff(data(i).velocity.time(1:2))));
    adcp.burst.cellsize = diff(data(i).velocity.y(1:2));
    adcp.burst.time = data(i).velocity.time;
    adcp.burst.range = data(i).velocity.y';
    adcp.burst.vel = cat(3,data(i).velocity.u,data(i).velocity.v,data(i).velocity.w);
    ax_vel = adcpQuickPlot(gcf,adcp,'vel',.4,NaT,[0 1],8);
    ylim(ax_vel(1),[0 max(adcp.burst.range)])

    figure(10*i+6); clf
    subplot(2,1,1)
    plot(data(i).velocity.time,data(i).velocity.pressure)
    subplot(2,1,2); hold on
    plot(data(i).velocity.time,data(i).velocity.pitch)
    % plot(data(i).velocity.time,data(i).velocity.roll)

end