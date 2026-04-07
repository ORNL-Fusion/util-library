function CoilGeometry = define_PROTO_coils
% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil 
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to z=0. The coil extends from        
% this location in the direction of larger z
%
% layers in radius, turns in axial
% coils centered around r=0
%
% Windings are constructed downstream in build_circular_coil*. The
% numbering used there is shown schematically below
%  r2 -------------------
%     |     |     |     |
%     |  2  |  4  |  6  |        % two layers
%     |     |     |     |
%     -------------------
%     |     |     |     |
%     |  1  |  3  |  5  |        % three turns
%     |     |     |     |
%  r1 -------------------
%    z1                  z1+dz
%
% JD Lore

CoilGeometry.ncoils = 13;
omat = ones(1,CoilGeometry.ncoils);
CoilGeometry.nturns  = 8*omat;
CoilGeometry.nlayers = 5*omat;
CoilGeometry.rr1 = 0.1221*omat;
CoilGeometry.rr2 = 0.1785*omat;
CoilGeometry.cl =  0.0979*omat;
CoilGeometry.z0 = [0.9392, 1.2493, 1.5792, 1.8151, 2.1413, 2.3392, 2.8953, 3.0603, 3.3398, 3.54289, 3.7335, 3.9240, 4.2414].';
