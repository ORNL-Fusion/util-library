function [Atheta] = afield_circular_coilset(coil,current,r,z)

nwind = length(coil.rwind);
Atheta = 0;
for i=1:nwind
    [Atheta0] = afield_circular_coil_analytic(coil.rwind(i),coil.zwind(i),r,z);
    Atheta = Atheta + Atheta0*current(i);
end
