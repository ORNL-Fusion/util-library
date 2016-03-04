function [df,ierr] = fl_derivs_dz_just_coils(Z,RP,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end

N = length(RP);
[Br,Bphi,Bz]=bfield_bs_cyl(RP(1:2:N-1),RP(2:2:N),Z*ones(size(RP(1:2:N-1))),bfield.coil,bfield.current,nowarn);
df(1:2:N-1) = Br./Bz; % dR/dz
df(2:2:N)   = Bphi./(Bz.*RP(1:2:N-1).'); % dphi/dz
ierr = 0;