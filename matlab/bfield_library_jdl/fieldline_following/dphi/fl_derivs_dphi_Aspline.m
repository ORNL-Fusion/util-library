function [df,ierr] = fl_derivs_dphi_Aspline(phi,RZ,bfield,nowarn)

N = length(RZ);
[Br,Bz,Bphi]=bfield_Aspline(RZ(1:2:N-1),RZ(2:2:N),phi*ones(size(RZ(1:2:N-1))),bfield.Arcoeff,bfield.Azcoeff,bfield.Aphicoeff,bfield.spline_info);              

df(1:2:N-1) = RZ(1:2:N-1).'.*Br./Bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bz./Bphi;
ierr = 0;