function [df,ierr] = fl_derivs_dl_MPEX(RPZ,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end
N = length(RPZ);

[Br,Bz] = bfield_circular_coil(bfield.coil,bfield.current,RPZ(1:3:N-2),RPZ(3:3:N));
Bphi = zeros(size(Br));
Btot = sqrt(Br.^2 + Bz.^2);

df(1:3:N-2) = Br./Btot; % dR/dl
df(2:3:N-1) = Bphi./(RPZ(1:3:N-2).'.*Btot); % dphi/dl
df(3:3:N)   = Bz./Btot; % dZ/dl
ierr = 0;