function [vel,qa] = interpolateADCP(vel,qa,t_max,r_max)
% Interpolate ADCP data removed by correlation filter.
%
%

[nt,nc,nb] = size(vel);

% loop through beams
for i = 1:nb
    veli = vel(:,:,i);

    % create duplicate with bad data replaced by nans
    veli_nan = veli;
    veli_nan(qa(:,:,i)) = nan;

    % interpolate cell by cell in time
    for j = 1:nc
        vj = veli(:,j);
        vj_nan = veli_nan(:,j);
        [blocks,widths] = findBlocks(isnan(vj_nan)); % find blocks of nans
        if isempty(blocks)
            continue
        end
        if blocks(1,1)==1 % get rid of blocks adjacent to beginning and end (cannot interpolate)
            blocks(1,:) = [];
            widths(1) = [];
        end
        if blocks(end,2)==nt
            blocks(end,:) = [];
            widths(end) = [];
        end

        % loop through blocks, interpolate if appropriate
        for k = 1:length(widths)
            if widths(k)<t_max % check that block is short enough to interpolate
                xk = (blocks(k,1)-1):(blocks(k,2)+1);
                vj(xk(2:end-1)) = interp1(xk([1 end]),vj_nan(xk([1 end])),xk(2:end-1));
                qa(xk(2:end-1),j,i) = false;
            end
        end
        veli(:,j) = vj;
    end
    vel(:,:,i) = veli;
end



