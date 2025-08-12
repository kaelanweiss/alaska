fid = fopen('F:/alaska2022/data/iceberg_surveys/raw/20220822_oyster/dragon/QGroundControl/2022-08-22 20-08-25 vehicle1.csv','r');
fgetl(fid);
line = fgetl(fid);
fclose(fid);
n = 50000;

tic;
for i = 1:n
    vals = strsplit(line,',');
    ts = vals{1};
end
t1 = toc;

tic;
for i = 1:n
    ts = line(1:find(line==',',1)-1);
end
t2 = toc;
