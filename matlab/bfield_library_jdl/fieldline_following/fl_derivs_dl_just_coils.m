function [df,ierr] = fl_derivs_dl_just_coils(RPZ,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end
N = length(RPZ);

[Br,Bphi,Bz]=bfield_bs_cyl(RPZ(1:3:N-2),RPZ(2:3:N-1),RPZ(3:3:N),bfield.coil,bfield.current,nowarn);
Btot = sqrt(Br.^2 + Bz.^2 + Bphi.^2);

df(1:3:N-2) = Br./Btot; % dR/dl
df(2:3:N-1) = Bphi./(RPZ(1:3:N-2).'.*Btot); % dphi/dl
df(3:3:N)   = Bz./Btot; % dZ/dl
ierr = 0;