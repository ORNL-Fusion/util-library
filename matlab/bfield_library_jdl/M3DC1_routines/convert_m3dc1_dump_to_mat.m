% convert_m3dc1_dump_to_mat
clearvars;

fname_in = 'C:\Work\M3DC1\Updated_148712_files_10_25_2013\m3dc1_148712_4101_1mm.dump';
fname_out = 'C:\Work\M3DC1\Updated_148712_files_10_25_2013\m3dc1_148712_4101_1mm.mat';


% % % fname_in = 'C:/Work/M3DC1/m3dc1_148712_4101.dump';
% % % fname_out = 'C:/Work/M3DC1/m3dc1_148712_4101.mat';
% % % fname_in = 'm3dc1_148712_4101_vac.dump';
% % % fname_out = 'm3dc1_148712_4101_vac.mat';

fid1 = fopen([fname_in,'1'],'r');
dat_line = fgetl(fid1);b.grid_path = sscanf(dat_line,'%s');
dat_line = fgetl(fid1);b.source = sscanf(dat_line,'%s');
dat_line = fgetl(fid1); data = sscanf(dat_line,'%i %i %i');
b.nr = data(1);b.nz = data(2);b.nn = data(3);
dat_line = fgetl(fid1); data = sscanf(dat_line,'%f %f %f');
b.factor = data(1); b.dr = data(2); b.dz = data(3);
b.npts = b.nr*b.nz;

%
% M3DC1 Fields
%
d2 = load([fname_in,'2']);

istart = 1; iend = b.nr;
b.r = d2(istart:iend);
istart = iend + 1; iend = istart + b.nz - 1;
b.z = d2(istart:iend);

istart = iend + 1; iend = istart + b.npts - 1;
b.br_amp = d2(istart:iend);
istart = iend + 1; iend = istart + b.npts - 1;
b.bt_amp = d2(istart:iend);
istart = iend + 1; iend = istart + b.npts - 1;
b.bz_amp = d2(istart:iend);

istart = iend + 1; iend = istart + b.npts - 1;
b.br_phase = d2(istart:iend);
istart = iend + 1; iend = istart + b.npts - 1;
b.bt_phase = d2(istart:iend);
istart = iend + 1; iend = istart + b.npts - 1;
b.bz_phase = d2(istart:iend);


%
% EQ Fields
%
d3 = load([fname_in,'3']);
istart = 1; iend = b.nr;
% b.r = d3(istart:iend);
istart = iend + 1; iend = istart + b.nz - 1;
% b.z = d3(istart:iend);

istart = iend + 1; iend = istart + b.npts - 1;
b.br_eq = d3(istart:iend);
istart = iend + 1; iend = istart + b.npts - 1;
b.bt_eq = d3(istart:iend);
istart = iend + 1; iend = istart + b.npts - 1;
b.bz_eq = d3(istart:iend);
close all;

save(fname_out,'b');

