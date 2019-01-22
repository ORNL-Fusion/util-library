clearvars;

fdir = 'C:\Work_archive\RUN_ARCHIVE\FLARE_RUNS\make_grid\W7X\OP12a_43kA_grid\';
fname_in = fullfile(fdir,'vessel.dat');
fname_out = fullfile(fdir,'vessel_new.dat');


fid = fopen(fname_in,'r');
if fid == -1
    error(['Could not open file: ',fullfile(div_path,divfname)])
end

% Just read lines until we get three ints
max_tries = 10;
icount = 0;
data = [];
dlabel = fgetl(fid);
while isempty(data)    
    dat_line = fgetl(fid);
    data = sscanf(dat_line,'%d %d %d %f %f %f %d');
    icount = icount + 1;
    if icount > max_tries
        error('max tries exceeded in read_limiter_file')
    end
end

ntor = data(1);
np_div = data(2);
ntper_lim = data(3); 
rshift = data(4)/100;
zshift = data(5)/100;
for it = 1:ntor
    dat_line = fgetl(fid);
    phivals(it) = (pi/180)*sscanf(dat_line,'%f');
    for jp = 1:np_div
        dat_line = fgetl(fid);
        data = sscanf(dat_line,'%f %f %d');
        div_r(it,jp) = data(1)/100 + rshift;
        div_z(it,jp) = data(2)/100 + zshift;        
    end
end
fclose(fid);


% Find 36 deg
[dphi,iphi] = min(abs(phivals-36*pi/180));



jp = 43:45;


figure; hold on; box on;
plot(div_r(iphi,:),div_z(iphi,:),'o-')
plot(div_r(iphi,jp),div_z(iphi,jp),'x-')

div_r(iphi,43) = 4.45; div_z(iphi,43) = -0.64;
div_r(iphi,44) = 4.5; div_z(iphi,44) = -0.7;
div_r(iphi,45) = 4.6; div_z(iphi,45) = -0.73;
plot(div_r(iphi,jp),div_z(iphi,jp),'x-')
plot(div_r(iphi,jp),-div_z(iphi,jp),'x-')

jp2 = 29:31;
plot(div_r(iphi,jp2),div_z(iphi,jp2),'x-')
div_r(iphi,jp2) = div_r(iphi,fliplr(jp));
div_z(iphi,jp2) = -div_z(iphi,fliplr(jp));
plot(div_r(iphi,:),div_z(iphi,:),'<-')
% asdfasdfds


% WRITE FILE



% dlabel = '# blah blah';
rshift = 0;  % in cm
zshift = 0;

fid = fopen(fname_out,'w');
if fid == -1
    error('Could not open file: %s\n',fname_out)
end

% write label and geom info
fprintf(fid,'%s\n',dlabel);
fprintf(fid,'%i %i %i %18.12f %18.12f\n', [ntor,np_div,ntper_lim,rshift,zshift]);

for it = 1:ntor
    fprintf(fid,'%18.12f\n',180*phivals(it)/pi);  % Phi in degrees
    for jp = 1:np_div
        fprintf(fid,'%18.12f %18.12f\n',100*div_r(it,jp),100*div_z(it,jp));  % toroidally symmetric R,Z in (cm)
    end
end

fclose(fid);