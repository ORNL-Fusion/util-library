function [coil,current] = build_Proto_coils_jackson(helicon_current,current_A,current_B,config,verbose,current_C)
if nargin < 4
    error('must specify inputs')
end


[nturns,nlayers,rr1,rr2,cl,z0,cur] = setup_Proto_coils(helicon_current,current_A,current_B,config,verbose,current_C);

ncoil = length(nturns);
for i = 1:ncoil
    coil_an = build_circular_coil_jackson(rr1(i),rr2(i),z0(i),cl(i),nturns(i),nlayers(i));
    if i > 1
        coil.rwind = [coil.rwind,coil_an.rwind];
        coil.zwind = [coil.zwind,coil_an.zwind];
        current = [current,cur(i)*ones(size(coil_an.rwind))];
    else
        coil.rwind = coil_an.rwind;
        coil.zwind = coil_an.zwind;
        current = cur(i)*ones(size(coil_an.rwind));
    end
end

