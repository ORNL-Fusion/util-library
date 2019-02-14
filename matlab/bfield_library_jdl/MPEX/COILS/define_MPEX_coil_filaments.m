function fil = define_MPEX_coil_filaments
% nturns = number of turns per coil layer in each coil (see note below)
% nlayers = number of layers in each coil (see note below)
% rr1 = inner radius of each coil
% rr2 = outer radius of each coil
% cl = axial length of each coil
% z0 = (axial) starting position of each coil relative to z=0 (here ECH center). The coil extends from        
% this location in the direction of larger z
%
% layers in radius, turns in axial
% coils centered around r=0
% windings ordered from min Z, min R; R increases fastest
%  r2 -------------------
%     |     |     |     |
%     |  4  |  5  |  6  |        % two layers
%     |     |     |     |
%     -------------------
%     |     |     |     |
%     |  1  |  2  |  3  |        % three turns
%     |     |     |     |
%  r1 -------------------
%    z1                  z1+dz

fil.ncoils = 21;
omat = ones(1,fil.ncoils);
fil.nturns  = 1*omat;
fil.nlayers = 1*omat;
fil.nwind = fil.nturns.*fil.nlayers;


% Coils ordered as [Helicon, ECH, ICH, Transport, Target]
% Number of coils in each: 6,4,5,3,3
fil.rr1 = [32.37*ones(1,6),32.37*ones(1,4),32.37*ones(1,5),67*ones(1,3)  ,67*ones(1,3)  ]*0.0254/2;  % Values are from diameter table
fil.thick = [1.5*ones(1,6),  1.5*ones(1,4),  1.5*ones(1,5),1.25*ones(1,3),1.25*ones(1,3)]*0.0254;
fil.rr2 = fil.rr1 + fil.thick;
fil.cl = [6*ones(1,6),6*ones(1,4),6*ones(1,5),15*ones(1,3),15*ones(1,3)]*0.0254;
fil.z0  = [-170, -158, -118, -104, -64, -52, -28, -14, 8, 19, 35, 49, 62.5, 76.5, 93.5, 111, 131, 151, 182, 202, 222]*0.0254;


fil.area = fil.cl.*fil.thick;