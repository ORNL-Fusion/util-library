% function test_circular_coil_field
clearvars;


[nturns,nlayers,rr1,rr2,cl,z0] = define_proto_coil_filaments;
% [nturns,nlayers,rr1,rr2,cl,z0,cur] = setup_Proto_coils;

% nturns = 1; nlayers = 1;
zmax = 5;    % axial distance to calculate the fields over
rmax = 0.4;  % radial distance to calculate the fields over
nr = 3;    % number of radial positions at which to calculate the field quantities
nz = 5;    % number of axial positions at which to calculate the field quantities

cur(:) = 0; cur(1) = 3300;

% [brg,bzg,atg,avec,zvec]=geom_coaxial_coilsm(z0,cl,rr1,rr2,nturns,nlayers,zmax,rmax,nz,nr);
% [br,bz,psi]=B_coaxial_coilsm(brg,bzg,atg,avec,cur);


% BIOT-SAVART

r1 = 0.1221;
r2 = 0.1785;
z1 = 0.9392;
dz = 0.0979;
nturns = 8; nlayers = 5;
% nturns = 1; nlayers = 1;
cur1 = 3300;
[coil,current] = build_circular_coil(r1,r2,z1,dz,nturns,nlayers,cur1);

% P_r = 1e-6;
% P_phi = 0;
% P_z = 1.15;
ir_test = 1
iz_test = 1
% P_r = avec(ir_test);
P_phi = 0;
% P_z = zvec(iz_test);
P_r = linspace(0.01,0.11,100);
P_z = (z1+dz/2)*ones(size(P_r));
[Br,Bphi,Bz]=bfield_bs_cyl(P_r,P_phi,P_z,coil,current);


% ANALYTICAL
coil_an = build_circular_coil_jackson(r1,r2,z1,dz,nturns,nlayers);
% [Br_an,Bz_an,Atheta_an] = bfield_circular_coil(coil_an,cur(1),P_r,P_z);

fprintf('-----------------------------------------\n')
fprintf('Evaluating at [R,Z] = %f,%f\n',P_r,P_z)
fprintf('Rick [br,bz] = %e,%e\n',br(ir_test,iz_test),bz(ir_test,iz_test))
fprintf('Biot [br,bz] = %e,%e\n',Br,Bz)
fprintf('Ban  [br,bz] = %e,%e\n',Br_an,Bz_an)
fprintf('-----------------------------------------\n')