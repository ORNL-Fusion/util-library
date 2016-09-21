function [df,ierr] = fl_derivs_dphi_xpand(phi,RZ,bfield,nowarn,field_choice)
% field_choice == 0: pert
% field_choice == 1: vacuum
N = length(RZ);
if field_choice == 1
    bfield.xpand.g = bfield.g;
end    

[Br,Bz,Bphi]=bfield_xpand(RZ(1:2:N-1),RZ(2:2:N),phi*ones(size(RZ(1:2:N-1))),bfield.xpand,nowarn,field_choice);

df(1:2:N-1) = RZ(1:2:N-1).'.*Br./Bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bz./Bphi;
ierr = 0;