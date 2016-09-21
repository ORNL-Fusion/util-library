clearvars;

run_path = 'C:\Work\DIII-D\164723\VMEC_XPAND\3059\';
fname = fullfile(run_path,'xpand_164723_3059.dat');

gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
g = readg_g3d(gfile_name);

xpand = read_xpand_field_file(fname);


TEST = 2;

phistart = 0;
Rstart = 2.2;
Zstart = 0.05;

dphi = 0.1*pi/180;
ntransits = 1;
nsteps = ntransits*2*pi/dphi;

if TEST == 1
    bfield.type = 'xpand_pert'; 
elseif TEST == 2
    bfield.type = 'xpand_vac'; 
    bfield.g =  g;
end
bfield.xpand = xpand;

f = follow_fieldlines_rzphi_dphi(bfield,Rstart,Zstart,phistart,dphi,nsteps);
x = f.r.*cos(f.phi);
y = f.r.*sin(f.phi);
z = f.z;


figure; hold on; box on;
plot3(x,y,z,'k')
