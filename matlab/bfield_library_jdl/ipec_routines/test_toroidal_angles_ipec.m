
clearvars;


run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\high\ipec\';
ipec = open_ipec_field(run_path);



neval = 100;
phi_eval = linspace(0,2*pi,neval);
Reval = 2.2*ones(size(phi_eval));
Zeval = 0.05*ones(size(phi_eval));
ifield_eval = 3;
[Br,Bz,Bphi,Btot]=bfield_ipec(Reval,Zeval,phi_eval,ipec,1,ifield_eval);
 rmp = build_d3d_icoils_jl([-2903.   2939.  -2889.   2935.  -2886.   2940. -2851.   2907.  -2866.   2918.  -2910.   2918.]);
[Br2,Bphi2,Bz2]=bfield_bs_cyl(Reval,phi_eval,Zeval,rmp.coil,rmp.current);
Btot2 = sqrt(Br2.^2 + Bphi2.^2 + Bz2.^2);


figure; hold on; box on;
plot(phi_eval,Br,'k','linewidth',2)
plot(phi_eval,Bz,'k--','linewidth',2)
plot(phi_eval,Bphi,'k-.','linewidth',2)
plot(phi_eval,Br2,'r','linewidth',2)
plot(phi_eval,Bz2,'r--','linewidth',2)
plot(phi_eval,Bphi2,'r-.','linewidth',2)
legend('IPEC Br','IPEC Bz','IPEC Bphi','B-S Br','B-S Bz','B-S Bphi')


asfdsafasf

phistart = 0;
Rstart = 2.2;
Zstart = 0.05;

dphi = 0.1*pi/180;
ntransits = 1;
nsteps = ntransits*2*pi/dphi;


% bfield.type = 'ipec_eq';
TEST = 3; bfield.type = 'ipec_vac'; ifield_eval = 3;  %%% THIS ONE IS FOLLOW PERTURBED LINE AND EVAL Bpert only!
% bfield.type = 'ipec_vac_only'; 

bfield.ipec = ipec;

f = follow_fieldlines_rzphi_dphi(bfield,Rstart,Zstart,phistart,dphi,nsteps);
x = f.r.*cos(f.phi);
y = f.r.*sin(f.phi);
z = f.z;
[Br,Bz,Bphi,Btot]=bfield_ipec(f.r,f.z,f.phi,ipec,1,ifield_eval);


figure(2); hold on; box on;
plot(f.phi,Btot,'k','linewidth',2)
plot(f.phi,Br,'k--','linewidth',2)
plot(f.phi,Bz,'k.-','linewidth',2)

% asfasdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gfile_name = 'C:\Work\DIII-D\164723\g164723.03059_410';
g = readg_g3d(gfile_name);
bfield2.g = g;
nowarn = 1;

if TEST == 1
    bfield2.type = 'gfile';
elseif TEST == 2 || TEST == 3
   rmp = build_d3d_icoils_jl([-2903.   2939.  -2889.   2935.  -2886.   2940. -2851.   2907.  -2866.   2918.  -2910.   2918.]); 
   bfield2.type = 'gfile+coils';
   bfield2.coil = rmp.coil;
   bfield2.current = rmp.current;
end


f2 = follow_fieldlines_rzphi_dphi(bfield2,Rstart,Zstart,phistart,dphi,nsteps);
x2 = f2.r.*cos(f.phi);
y2 = f2.r.*sin(f.phi);
z2 = f2.z;

[Bout,ierr] = bfield_geq_bicub(g,f2.r,f2.z,nowarn);
if TEST == 2 
    [Br,Bphi,Bz]=bfield_bs_cyl(f2.r,f2.phi,f2.z,bfield2.coil,bfield2.current,nowarn);
    Bout.br = Bout.br + Br;
    Bout.bphi = Bout.bphi + Bphi;
    Bout.bz = Bout.bz + Bz;
elseif TEST == 3
    [Br,Bphi,Bz]=bfield_bs_cyl(f2.r,f2.phi,f2.z,bfield2.coil,bfield2.current,nowarn);
    Bout.br = Br;
    Bout.bphi = Bphi;
    Bout.bz = Bz;    
end

Btot_2 = sqrt(Bout.br.^2 + Bout.bphi.^2 + Bout.bz.^2);

figure(2);
plot(f2.phi,Btot_2,'r')
plot(f2.phi,Bout.br,'r--')
plot(f2.phi,Bout.bz,'r.-')