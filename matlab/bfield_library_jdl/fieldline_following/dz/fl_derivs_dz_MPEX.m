function [df,ierr] = fl_derivs_dz_MPEX(Z,RP,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end

N = length(RP);
[Br,Bz] = bfield_circular_coils(bfield.coil,bfield.current,RP(1:2:N-1),Z);
df(1:2:N-1) = Br./Bz; % dR/dz
df(2:2:N)   = zeros(size(Br)); % dphi/dz
ierr = 0;