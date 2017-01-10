function [df,ierr] = fl_derivs_dphi_Bgrid(phi,RZ,bfield,nowarn)

N = length(RZ);
[Br,Bz,Bphi]=bfield_grid(RZ(1:2:N-1),RZ(2:2:N),phi*ones(size(RZ(1:2:N-1))),bfield.Bgrid,nowarn);

df(1:2:N-1) = RZ(1:2:N-1).'.*Br./Bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bz./Bphi;
ierr = 0;