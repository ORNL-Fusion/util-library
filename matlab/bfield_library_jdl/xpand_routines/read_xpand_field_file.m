function xpand = read_xpand_field_file(fname)
% xpand.r(ir,iz) -- indices start from lower left corner

% clearvars;
% run_path = 'C:\Work\DIII-D\164723\VMEC_XPAND\3059\';
% fname = fullfile(run_path,'xpand_164723_3059.dat');

[base,fname_part,ext] = fileparts(fname);
matfile = fullfile(base,[fname_part,ext,'.mat']);
if exist(matfile,'file')
    fprintf('Reading xpand field data from .mat file\n')
    load(matfile);
    return;
else
    fprintf('Reading xpand data from ASCII file\n')
end

% Read header info


fid = fopen(fname);
dat = fgetl(fid);
fprintf('%s\n',dat)
dat = fgetl(fid);  % NR = 
dat2 = sscanf(dat,'%s %s %s %d');
nr = dat2(5);
fprintf('NR = %d\n',nr)

dat = fgetl(fid);  % Np = 
dat2 = sscanf(dat,'%s %s %s %d');
nphi = dat2(5);
fprintf('Nphi = %d\n',nphi)

dat = fgetl(fid);  % NZ = 
dat2 = sscanf(dat,'%s %s %s %d');
nz = dat2(5);
fprintf('NZ = %d\n',nz)

dat = fgetl(fid); %label
fclose(fid);

% Read bulk data
data = dlmread(fname,'',5,0);

% Reshape array
xpand.r = data(1:nr,1);
xpand.z = data(1:nr:nr*nz,3);
xpand.phi = data(1:nr*nz:nphi*nr*nz,2);

xpand.nphi = nphi;
xpand.nz = nz;
xpand.nr = nr;

% r_tmp = data(:,1);
% rgood = reshape(r_tmp,[nr,nz,nphi]);
% z_tmp = data(:,3);
% zgood = reshape(z_tmp,[nr,nz,nphi]);
% phi_tmp = data(:,2);
% phigood = reshape(phi_tmp,[nr,nz,nphi]);


xpand.Br   = reshape(data(:,4),[nr,nz,nphi]);
xpand.Bphi = reshape(data(:,5),[nr,nz,nphi]);
xpand.Bz   = reshape(data(:,6),[nr,nz,nphi]);
xpand.pres    = reshape(data(:,7),[nr,nz,nphi]);
xpand.Brvac   = reshape(data(:,8),[nr,nz,nphi]);
xpand.Bphivac = reshape(data(:,9),[nr,nz,nphi]);
xpand.Bzvac   = reshape(data(:,10),[nr,nz,nphi]);

% Add first slice to end to get periodic grid
xpand.nphi = xpand.nphi + 1;
xpand.Br(:,:,xpand.nphi) = xpand.Br(:,:,1);
xpand.Bphi(:,:,xpand.nphi) = xpand.Bphi(:,:,1);
xpand.Bz(:,:,xpand.nphi) = xpand.Bz(:,:,1);
xpand.pres(:,:,xpand.nphi) = xpand.pres(:,:,1);
xpand.Brvac(:,:,xpand.nphi) = xpand.Brvac(:,:,1);
xpand.Bphivac(:,:,xpand.nphi) = xpand.Bphivac(:,:,1);
xpand.Bzvac(:,:,xpand.nphi) = xpand.Bzvac(:,:,1);

xpand.phi(xpand.nphi) = 2*pi;

save(matfile,'xpand')
