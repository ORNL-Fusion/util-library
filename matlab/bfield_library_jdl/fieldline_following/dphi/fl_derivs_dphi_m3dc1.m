function [df,ierr] = fl_derivs_dphi_m3dc1(phi,RZ,rmp,nowarn)
if nargin < 4
    nowarn = 0;
end
N = length(RZ);

[Br,Bz,Bphi,ierr] = bfield_m3dc1(rmp.b,RZ(1:2:N-1),phi,RZ(2:2:N),rmp.scale);
Bout.br = Br;
Bout.bphi = Bphi;
Bout.bz = Bz;

if ierr == 1
    error('Error in calling bfield_m3dc1!')
end

df(1:2:N-1) = RZ(1:2:N-1).'.*Bout.br./Bout.bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bout.bz./Bout.bphi;
ierr = 0;