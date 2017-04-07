function fil = define_proto_coil_filaments
% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil (see note below)
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to the start of the magnetic field profile at z=0. The coil extends from        
% this location in the direction of larger z
fil.ncoils = 12;
omat = ones(1,fil.ncoils);
fil.nturns  = 8*omat;
fil.nlayers = 5*omat;
fil.rr1 = 0.1221*omat;
fil.rr2 = 0.1785*omat;
fil.cl = 0.0979*omat;
fil.z0  = [0.9392    1.2492    1.5792    1.8152    2.1412    2.3392    2.8952    3.1712    3.3692    3.6852    3.9992    4.3172];


