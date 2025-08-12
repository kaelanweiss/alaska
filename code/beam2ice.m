function burst = beam2ice(burst,beam_distance)
% Function to map distance-from-ADCP data to distance-from-ice data.
%
% burst = beam2ice(burst,beam_distance)
%
% KJW
% 27 Sep 2022

FLDS = {'vel','cor','amp'};
nb = length(beam_distance);
nt = length(burst.time);

% find furthest extent, calculate r(t)
r_all = [];
for i = 1:nb
    r_all = [r_all; beam_distance(i).r];
    beam_distance(i).r_t = interp1(beam_distance(i).t,beam_distance(i).r,burst.time,'previous','extrap')+burst.cellsize;
end

% count how many cells we'll need, build distance-from-ice range
nc = ceil(max(r_all-burst.blankingdist)/burst.cellsize)+1;
r_ice = ((1:nc)-1)*burst.cellsize;

% loop through fields
for i = 1:length(FLDS)
    fldi = FLDS{i};
    tmp = NaN(nt,nc,nb); % preallocate
    % loop through beams
    for j = 1:nb
        % loop through time
        for k = 1:nt
            % find slice in original, flip it, put it in output (tmp)
            r_idx = burst.range <= beam_distance(j).r_t(k);
            prof = flip(squeeze(burst.(fldi)(k,r_idx,j)));
            tmp(k,1:length(prof),j) = prof;
        end
    end
    burst.(fldi) = tmp;
end
burst.range = r_ice;
        