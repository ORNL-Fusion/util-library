function [nturns,nlayers,rr1,rr2,cl,z0] = define_proto_coil_filaments
% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil (see note below)
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to the start of the magnetic field profile at z=0. The coil extends from        this location in the direction of larger z
ncoils = 12;
omat = ones(1,ncoils);
nturns  = 8*omat;
nlayers = 5*omat;
rr1 = 0.1221*omat;
rr2 = 0.1785*omat;
cl = 0.0979*omat;
z0  = [0.9392    1.2492    1.5792    1.8152    2.1412    2.3392    2.8952    3.1712    3.3692    3.6852    3.9992    4.3172];
