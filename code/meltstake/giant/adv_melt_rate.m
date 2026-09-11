load F:\GIANT2026\meltstake\data\proc\20260803_1430_ms01\adv\pck_edges.mat

nt = length(edges.time);

edges.time_hr = hours(edges.time-edges.time(1));

