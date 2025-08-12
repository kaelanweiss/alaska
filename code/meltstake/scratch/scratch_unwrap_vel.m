% run "ms_vel_scale_osm.m" with dep_num 16 section 2 to retrieve data

idxt = true(size(adcp.burst.time));
time = adcp.burst.time;

% get raw beam velocity and correlation
beam_num = 2;
cell_num = 1;
v = adcp.burst.vel(idxt,cell_num,beam_num);
c = adcp.burst.cor(idxt,cell_num,beam_num);

% empirical velocity range
vel_range = diff(extrema(v));

% qa a lil bit
v(c<30) = nan;

% diff in time
dv = diff(v);

% find jumps
jump_threshold = 0.2;
jump = fix(dv/jump_threshold);
idxj = jump ~= 0 & ~isnan(jump);
jump_inds = find(idxj);

% plot some stuff
clear ax
figure(1); clf

ax(1) = subplot(4,1,1:2);
hold on
plot(extrema(time),max(abs(v))*[1 1],'k--')
plot(extrema(time),-max(abs(v))*[1 1],'k--')
plot(time,v,'.-','color',colors(1))
plot(time(idxj),v(idxj),'ro')
ylim(.5*[-1 1])

ax(2) = subplot(4,1,3);
plot(time,[dv; 0],'')

ax(3) = subplot(4,1,4);
plot(time,[jump; 0])

linkaxes(ax,'x')

% start fixin stuff
v_fix1 = v;
for i = 1:length(jump_inds)-1
    j1 = jump(jump_inds(i));
    j2 = jump(jump_inds(i+1));

    if j1*j2 < 0
        slice = jump_inds(i)+1:jump_inds(i+1);
        if any(isnan(v(slice)))
            continue
        end
        if j1 < 0
            v_fix1(slice) = v_fix1(slice) + vel_range;
        else
            v_fix1(slice) = v_fix1(slice) - vel_range;
        end
    end
end

plot(ax(1),time,v_fix1+.01,'.-','color',colors(2))

% let's try this a different way
tic
v_fix2 = v;
nt = length(v);
i = 1;
count = 0;

while i < nt
    count = count + 1;
    v_rem = v_fix2(i:end);
    jump_rem = fix(diff(v_rem)/jump_threshold);
    j_inds = find(abs(jump_rem)>jump_threshold,2);

    if isempty(j_inds)
        break
    end

    j1 = jump_rem(j_inds(1));
    j2 = jump_rem(j_inds(2));

    if j1*j2 < 0
        slice = j_inds(1)+1:j_inds(2);
        if any(isnan(v_rem(slice)))
            i = j_inds(2) + i - 1;
            continue
        end
        
        % put fixed slice back in original
        if j1 < 0
            v_fix2(slice+i-1) = v_rem(slice) + vel_range;
        else
            v_fix2(slice+i-1) = v_rem(slice) - vel_range;
        end
    end
    
    i = j_inds(2) + i - 1;
    fprintf('%d/%d\n',i,nt)
end
toc

plot(ax(1),time,v_fix2+.02,'.-','color',colors(5))