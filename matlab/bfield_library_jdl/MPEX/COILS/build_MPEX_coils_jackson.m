function [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson(current_in,verbose)
%TODO: Add description
% J.D. Lore
if nargin == 0
    error('Must specify inputs')
end
if nargin < 2
    verbose = 0;
end

[CoilGeometry,currentPerWinding] = setup_MPEX_coils(current_in,verbose);
nwind = CoilGeometry.nturns.*CoilGeometry.nlayers;
nwind_tot = sum(CoilGeometry.nturns.*CoilGeometry.nlayers);
windingCurrent = zeros(1,nwind_tot);
Coil.rwind = zeros(1,nwind_tot);
Coil.zwind = zeros(1,nwind_tot);
for i = 1:CoilGeometry.ncoils
    coil_an = build_circular_coil_jackson(CoilGeometry.rr1(i),CoilGeometry.rr2(i),CoilGeometry.z0(i),CoilGeometry.cl(i),CoilGeometry.nturns(i),CoilGeometry.nlayers(i));
    i1 = 1 + sum(nwind(1:i-1));
    i2 = sum(nwind(1:i));
    Coil.rwind(i1:i2) = coil_an.rwind;
    Coil.zwind(i1:i2) = coil_an.zwind;
    windingCurrent(i1:i2) = currentPerWinding(i);
end
