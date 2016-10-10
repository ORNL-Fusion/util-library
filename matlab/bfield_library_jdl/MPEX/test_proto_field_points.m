% function test_proto_field_points
clearvars;

% [nturns,nlayers,rr1,rr2,cl,z0,cur] = setup_Proto_coils;
[nturns,nlayers,rr1,rr2,cl,z0] = define_proto_coil_filaments;

rmax = 0.1;
zmax = 5;
nz = 10;
nr = 5;
[brg,bzg,atg,avec,zvec]=geom_coaxial_coilsm(z0,cl,rr1,rr2,nturns,nlayers,zmax,rmax,nz,nr);

[br,bz,psi]=B_coaxial_coilsm(brg,bzg,atg,avec,cur);

% [br,bz,psi]=B_coaxial_coilsm(brg,bzg,atg,avec,cur);