function [Br,Bz,Atheta] = bfield_circular_coil(coil,current,r,z)

nwind = length(coil.rwind);

Br = 0;
Bz = 0;
Atheta = 0;
for i=1:nwind
    [Br0,Bz0,Atheta0] = bfield_MPEX_coils(coil.rwind(i),coil.zwind(i),r,z);
    Br = Br + Br0*current(i);
    Bz = Bz + Bz0*current(i);
    Atheta = Atheta + Atheta0*current(i);
end
