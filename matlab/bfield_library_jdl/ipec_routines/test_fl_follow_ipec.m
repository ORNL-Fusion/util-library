clearvars;

run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\low\gpec\';
ipec = open_ipec_field(run_path);

TEST = 2;

phistart = 0;
Rstart = 2.2;
Zstart = 0.05;

dphi = 0.1*pi/180;
ntransits = 1;
nsteps = ntransits*2*pi/dphi;

if TEST == 1
    bfield.type = 'ipec_eq'; 
elseif TEST == 2
    bfield.type = 'ipec_vac'; 
end
bfield.ipec = ipec;

f = follow_fieldlines_rzphi_dphi(bfield,Rstart,Zstart,phistart,dphi,nsteps);
x = f.r.*cos(f.phi);
y = f.r.*sin(f.phi);
z = f.z;


figure; hold on; box on;
plot3(x,y,z,'k')

% asfasdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear bfield
gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
g = readg_g3d(gfile_name);
bfield.g = g;

if TEST == 1
    bfield.type = 'gfile';
elseif TEST == 2
   rmp = build_d3d_icoils_jl([-2903.   2939.  -2889.   2935.  -2886.   2940. -2851.   2907.  -2866.   2918.  -2910.   2918.]); 
   bfield.type = 'gfile+coils';
   bfield.coil = rmp.coil;
   bfield.current = rmp.current;
end


f2 = follow_fieldlines_rzphi_dphi(bfield,Rstart,Zstart,phistart,dphi,nsteps);
x2 = f2.r.*cos(f.phi);
y2 = f2.r.*sin(f.phi);
z2 = f2.z;


% figure; hold on; box on;
plot3(x2,y2,z2,'r')