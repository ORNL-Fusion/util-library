clearvars;

fname = 'C:\Work\DIII-D\APS 2016\irtv01_156855_2000_4900.mat'; shot = 156855; elm_file = 'elm_156855_fs6middaf.mat';
load(fname);

thresh = 0.5;
dt_thresh = 10; %ms

twin = [4200,4800];
% twin = [3200,3800];
% elm_win = [0.8,0.95];
elm_win = [0.1,0.9];

DIV_PLOT = 1; % 1 = outer
if DIV_PLOT == 1
    drwant_min = -10;
    drwant_max = 20;
else
    error('set up inner divertor')
end

% RAW
time = double(ir.time);
data = double(ir.data);
drout = double(ir.drout);
drin = double(ir.drin);

% TIME SELECT
t_good = find(time >= twin(1) & time <= twin(2));
% eskeep = es(t_good);
time = time(t_good);
data = data(:,t_good);
drout = drout(:,t_good);
drin = drin(:,t_good);


drfind_elm = -10;

for i = 1:length(time)
    [dx,ix] = min(abs(drout(:,i) - drfind_elm));
    dat_1d(i) = data(ix,i);
end

sig = dat_1d;


spikes = find_spikes(time,sig,thresh,dt_thresh);



