% script to transform spider beam coordinate velocities to an
% ice-referenced earth coordinate system
%
% steps:
% 1) bin map velocity based on ice distance
% 2) transform to instrument coordinates
% 3) rotate to ice-referenced frame

clear

load F:alaska2022\data\iceberg_surveys\mat\20220824_singingflower\spider\adcp.mat
load F:/alaska2022/data/iceberg_surveys/mat/20220824_singingflower/spider/adcp_exclude.mat
load edge_idx_all

r = adcp.burst.range;
r_ice = r(edge_idx_all);
[nt,nc,nb] = size(adcp.burst.vel);

%% qa and exclusion
qa_cor = adcp.burst.cor < 60;
adcp.burst.vel(qa_cor) = nan;

exclude_idx(1:find(~exclude_idx,1,'first')) = 0;
adcp.burst.vel(repmat(exclude_idx,[1 nc nb])) = nan;

figure(1); clf
ax = adcpQuickPlot(figure(1),adcp,'vel',0.05*[-1 1],[0 Inf],[0 1],1);
for i = 1:4
    hold(ax(i),'on')
    plot(ax(i),adcp.burst.time,r_ice(:,i),'k')
end
%% bin map
%profile on
% find closest beam
[r_min,beam0] = min(r_ice,[],2);
nr = sum(r<=max(r_min))+1;

% set pivot for coordinate scaling
r0 = 0;

% squeeze each beam axis between the pivot and the ice distance
r_shift = 0*adcp.burst.vel;
for i = 1:nt
    ri = r_min(i);
    for j = 1:nb
        r_shift(i,:,j) = (r-r0)*(ri/r_ice(i,j)) + r0;
    end
end

% interpolate beam velocities onto closest beam scale
r_int = adcp.burst.range(1:nr);
vel_int = nan*adcp.burst.vel(:,1:nr,:);

for i = 1:nt
    % get range axis
    ri = r_shift(i,1:nr,beam0(i));
    % interpolate velocities
    %vel_int(i,:,beam0) = adcp.burst.vel(i,:,beam0);
    %j_list = 1:nb;
    %j_list(beam0(i)) = [];
    for j = 1:nb
        vj = interp1(r_shift(i,:,j),adcp.burst.vel(i,:,j),ri,'linear');
        vj(isnan(vj)) = interp1(r_shift(i,:,j),adcp.burst.vel(i,:,j),ri(isnan(vj)),'nearest');
        vel_int(i,:,j) = vj;
    end
end



%profile viewer

sum(isnan(vel_int),'all')

adcp2 = adcp;
adcp2.burst.vel = vel_int;
adcp2.burst.range = r_int;
adcpQuickPlot(figure(2),adcp2,'vel',0.05*[-1 1],[0 Inf],[0 1],1);

%% transform to instrument coordinates
% smooth with Hanning convolution
vel_smth = vel_int;
k = 8*60+1;
for i = 1:nr
    for j = 1:nb
        vel_smth(:,i,j) = hannFilter(vel_int(:,i,j),k);
    end
end

adcp3 = adcp2;
adcp3.burst.vel = vel_smth;
adcpQuickPlot(figure(3),adcp3,'vel',0.05*[-1 1],[0 Inf],[0 1],1);

% beam 2 instrument
vel_xyz = vel_smth;
rot180z = diag([1 -1 -1 -1]);
for i = 1:nt
    for j = 1:nr
        vel_xyz(i,j,:) = rot180z*adcp.burst.beam2xyz*squeeze(vel_smth(i,j,:));
    end
end

vz = (vel_xyz(:,:,3)+vel_xyz(:,:,4))/2;
ve = vel_xyz(:,:,3)-vel_xyz(:,:,4);

vel_xyz(:,:,3) = vz;
vel_xyz(:,:,4) = ve;

%%
adcp4 = adcp2;
adcp4.burst.vel = vel_xyz;
[ax,cbar] = adcpQuickPlot(figure(4),adcp4,'vel',0.05*[-1 1],[0 Inf],[0 Inf],10);
cbar_pos = cbar.Position;
delete(cbar)

lbls = {'right','up','away','error'};
clims = [6 6 2 2];
for i = 1:4
    ax(i).CLim = clims(i)*[-1 1]/100;
    cbar(i) = colorbar(ax(i));
    cbar(i).Position = [cbar_pos(1) ax(i).Position(2) cbar_pos(3) ax(i).Position(4)];
    cbar(i).Label.String = 'm/s';
    cbar(i).Label.FontSize = 11;
    text(ax(i),0.01,0.95,lbls{i},'units','normalized','fontsize',12)
end

