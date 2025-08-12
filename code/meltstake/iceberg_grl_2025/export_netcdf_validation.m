% Script to validate that the data in the exported netCDF files is correct
% and fully usable.

clear

% netcdf location, list files
nc_dir = 'F:meltstake/data/netcdf';
d = dir(fullfile(nc_dir,'XeitlSit*.nc'));
nf = length(d);

% summary output
fid = fopen(fullfile(nc_dir,'validation','segment_summary.txt'),'w');

% plot things
fs  = 11;
lw = 1;
figsize = [16 22];

% padding
pad = [.08 .05 .05 .1];
shift = [.02 0.05];

% loop through files
for i = 1:nf
    finfo = ncinfo(fullfile(d(i).folder,d(i).name));
    dimnames = {finfo.Dimensions.Name};
    allnames = {finfo.Variables.Name};
    varnames = setxor(allnames,dimnames);

    clear data dims
    data = struct();
    dims = struct();

    % load dims
    for j = 1:length(dimnames)
        posj = strcmp(allnames,dimnames{j});
        dims.(dimnames{j}) = struct();
        dims.(dimnames{j}).val = ncread(fullfile(d(i).folder,d(i).name),dimnames{j});
        if  ~isempty(finfo.Variables(posj).Attributes) && any(strcmp({finfo.Variables(posj).Attributes.Name},'units'))
            dims.(dimnames{j}).units = ncreadatt(fullfile(d(i).folder,d(i).name),dimnames{j},'units');
        else
            dims.(dimnames{j}).units = '';
        end
        % convert time dims to datetime
        if length(dimnames{j}) > 5 && strcmp(dimnames{j}(1:5),'time_')
            dims.(dimnames{j}).val = datetime(dims.(dimnames{j}).val,'convertFrom','posixtime');
        end
    end
    
    % load data
    for j = 1:length(varnames)
        posj = strcmp(allnames,varnames{j});
        data.(varnames{j}) = struct();
        data.(varnames{j}).dims = {finfo.Variables(posj).Dimensions.Name};
        data.(varnames{j}).val = ncread(fullfile(d(i).folder,d(i).name),varnames{j});
        if ~isempty(finfo.Variables(posj).Attributes) && any(strcmp({finfo.Variables(posj).Attributes.Name},'units'))
            data.(varnames{j}).units = ncreadatt(fullfile(d(i).folder,d(i).name),varnames{j},'units');
        else
            data.(varnames{j}).units = '';
        end
    end

    % print and plot information
    % print to file
    fprintf(fid,'Segment: %02d\n',dims.segment_number.val);
    fprintf(fid,'\tdepth: %.1f %s\n',round(dims.meltstake_depth.val,1),dims.meltstake_depth.units);
    fprintf(fid,'\tstart: %s\n',dims.start_time.val);
    fprintf(fid,'\tend: %s\n\n',dims.end_time.val);
    
    % plot
    fig = figure(1); clf
    setFigureSize(fig,figsize);
    clear ax

    % T
    ax(1) = axes(fig,'position',axgridpos(4,2,1,pad,shift));
    plot(dims.(data.T.dims{1}).val,data.T.val,'linewidth',lw)
    ylabel(data.T.units,'fontsize',fs)

    % S
    ax(2) = axes(fig,'position',axgridpos(4,2,2,pad,shift));
    hold(ax(2),'on')
    axis ij
    box on
    for j = 1:length(dims.time_S.val)
        plot(data.S.val(j,:),dims.z_S.val,'.-','linewidth',lw,'color',colors(2))
    end
    xlabel(data.S.units,'fontsize',fs)
    ylabel(dims.z_S.units,'fontsize',fs)

    % u
    j0 = 116;
    CLIM = max(abs([min([min(data.u.val,[],'all','omitnan') min(data.v.val,[],'all','omitnan') min(data.w.val,[],'all','omitnan')]),...
            max([max(data.u.val,[],'all','omitnan') max(data.v.val,[],'all','omitnan') max(data.w.val,[],'all','omitnan')])]));
    for j = 1:3
        ax(3+2*(j-1)) = axes(fig,'position',axgridpos(4,2,3+2*(j-1),pad,shift));
        pcolor(dims.time_uvw.val,dims.y_uvw.val,data.(char(j0+j)).val')
        shading flat
        cmocean('bal')
        clim(CLIM*[-1 1])
        ylabel(dims.y_uvw.units,'fontsize',fs)
        ylim([0 max(ylim)])
    end

    % echo
    CLIM = [min([min(data.echo1.val,[],'all','omitnan') min(data.echo2.val,[],'all','omitnan') min(data.echo3.val,[],'all','omitnan')]),...
            max([max(data.echo1.val,[],'all','omitnan') max(data.echo2.val,[],'all','omitnan') max(data.echo3.val,[],'all','omitnan')])];
    for j = 1:3
        ax(4+2*(j-1)) = axes(fig,'position',axgridpos(4,2,4+2*(j-1),pad,shift));
        pcolor(dims.time_echo.val,dims.dist_echo.val,data.(['echo' num2str(j)]).val')
        shading flat
        clim(CLIM)
        ylabel(dims.dist_echo.units,'fontsize',fs)
    end

    linkaxes(ax(setxor(1:8,2)),'x')
    xlim(ax(1),dims.time_T.val([1 end]))

    cbu = colorbar(ax(7),'orientation','horizontal','position',...
        [ax(7).Position(1) ax(7).Position(2)-0.08 ax(7).Position(3) 0.02]);
    cbu.Label.String = data.u.units;
    cbu.Label.FontSize = fs;

    cbe = colorbar(ax(8),'orientation','horizontal','position',...
        [ax(8).Position(1) ax(8).Position(2)-0.08 ax(8).Position(3) 0.02]);
    cbe.Label.String = data.echo1.units;
    cbe.Label.FontSize = fs;

    title(ax(1),sprintf('Segment %02d',dims.segment_number.val),'fontsize',fs+2)
    title(ax(2),dims.start_time.val,'fontsize',fs+2)
    
    remap_vars = [2 1 6 3 7 4 8 5];
    for j = 1:8
        txt = text(ax(j),1-.03,.97,varnames{remap_vars(j)},'units','normalized','fontsize',fs,...
            'verticalalignment','top','horizontalalignment','right');
    end

    % print figure
    print(figure(1),fullfile(nc_dir,'validation',sprintf('segment_%02d.png',dims.segment_number.val)),'-dpng','-r250')
end

fclose(fid);