function [df,ierr] = fl_derivs_dphi_just_coils(phi,RZ,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end
N = length(RZ);

[Br,Bphi,Bz]=bfield_bs_cyl(RZ(1:2:N-1),phi*ones(size(RZ(1:2:N-1))),RZ(2:2:N),bfield.coil,bfield.current,nowarn);

df(1:2:N-1) = RZ(1:2:N-1).'.*Br./Bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bz./Bphi;
ierr = 0;