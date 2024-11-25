function bfield = load_bgrid(bfieldDir,bfieldPrefix)

fnameMat = fullfile(bfieldDir,strcat(bfieldPrefix,'.mat'));
if exist(fnameMat)
    fprintf('Loading mat version %s\n',fnameMat)
    load(fnameMat);
else
    fprintf('Loading dat version %s\n',bfieldPrefix)
    Bgrid = loadBfieldFromDat(bfieldDir,bfieldPrefix);
end

bfield.type = 'Bgrid';
bfield.Bgrid = Bgrid;
bfield.nsym  = Bgrid.nsym;
end

function Bgrid = loadBfieldFromDat(bfieldDir,bfieldPrefix)
fName = fullfile(bfieldDir,strcat(bfieldPrefix,'_layout.dat'));

data = dlmread(fName);
Bgrid.nr = data(1);
Bgrid.nz = data(2);
Bgrid.nphi = data(3) + 1; % Will add in last plane below
Bgrid.nsym = data(4);
Bgrid.Rmin = data(5);
Bgrid.Rmax = data(6);
Bgrid.Zmin = data(7);
Bgrid.Zmax = data(8);
Bgrid.R = linspace(Bgrid.Rmin,Bgrid.Rmax,Bgrid.nr);
Bgrid.Z = linspace(Bgrid.Zmin,Bgrid.Zmax,Bgrid.nz);

Bgrid.dR = (Bgrid.Rmax - Bgrid.Rmin)/(Bgrid.nr-1);
Bgrid.dZ = (Bgrid.Zmax - Bgrid.Zmin)/(Bgrid.nz-1);
Bgrid.dphi = (2*pi/Bgrid.nsym)/(Bgrid.nphi - 1);
Bgrid.phi = linspace(0,2*pi/Bgrid.nsym,Bgrid.nphi);

fprintf('loading phi\n')
fName = fullfile(bfieldDir,strcat(bfieldPrefix,'_phi.dat'));
data = dlmread(fName);
Bgrid.Bphi = reshape(data, [Bgrid.nphi-1, Bgrid.nz, Bgrid.nr]);
Bgrid.Bphi = permute(Bgrid.Bphi, [3, 2, 1]);
Bgrid.Bphi(:,:,end+1) = Bgrid.Bphi(:,:,1);

fprintf('loading r\n')
fName = fullfile(bfieldDir,strcat(bfieldPrefix,'_r.dat'));
data = dlmread(fName);
Bgrid.Br = reshape(data, [Bgrid.nphi-1, Bgrid.nz, Bgrid.nr]);
Bgrid.Br = permute(Bgrid.Br, [3, 2, 1]);
Bgrid.Br(:,:,end+1) = Bgrid.Br(:,:,1);

fprintf('loading z\n')
fName = fullfile(bfieldDir,strcat(bfieldPrefix,'_z.dat'));
data = dlmread(fName);
Bgrid.Bz = reshape(data, [Bgrid.nphi-1, Bgrid.nz, Bgrid.nr]);
Bgrid.Bz = permute(Bgrid.Bz, [3, 2, 1]);
Bgrid.Bz(:,:,end+1) = Bgrid.Bz(:,:,1);

grid_file = fullfile(bfieldDir,strcat(bfieldPrefix,'.mat'));
save(grid_file,'Bgrid')
fprintf('saved file: %s\n',grid_file)

end