function [df,ierr] = fl_derivs_dphi_just_coils(phi,RZ,rmp,nowarn)
if nargin < 4
    nowarn = 0;
end
N = length(RZ);

[Br,Bphi,Bz]=bfield_bs_cyl(RZ(1:2:N-1),phi,RZ(2:2:N),rmp.coil,rmp.current,nowarn);
Bout.br = Br;
Bout.bphi = Bphi;
Bout.bz = Bz;

   
df(1:2:N-1) = RZ(1:2:N-1).'.*Bout.br./Bout.bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bout.bz./Bout.bphi;
ierr = 0;