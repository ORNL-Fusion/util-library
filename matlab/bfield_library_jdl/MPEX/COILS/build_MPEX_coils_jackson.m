function [coil,current] = build_MPEX_coils_jackson(current_in,verbose)
% J.D. Lore
if nargin == 0
    error('Must specify inputs')
end
if nargin < 2
    verbose = 0;
end
    
[fil,cur] = setup_MPEX_coils(current_in,verbose);
nwind = fil.nturns.*fil.nlayers;
nwind_tot = sum(fil.nturns.*fil.nlayers);
current = zeros(1,nwind_tot);
coil.rwind = zeros(1,nwind_tot);
coil.zwind = zeros(1,nwind_tot);
for i = 1:fil.ncoils
    coil_an = build_circular_coil_jackson(fil.rr1(i),fil.rr2(i),fil.z0(i),fil.cl(i),fil.nturns(i),fil.nlayers(i));
    i1 = 1 + sum(nwind(1:i-1));
    i2 = sum(nwind(1:i));
    coil.rwind(i1:i2) = coil_an.rwind;
    coil.zwind(i1:i2) = coil_an.zwind;
    current(i1:i2) = cur(i);
end
