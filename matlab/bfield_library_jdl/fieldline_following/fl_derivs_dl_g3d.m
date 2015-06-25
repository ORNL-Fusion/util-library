function [df,ierr] = fl_derivs_dl_g3d(L,RPZ,g,rmp,nowarn)
if nargin < 5
    nowarn = 0;
end
N = length(RPZ);

% Axisym part
[Bout,ierr_Bas] = bfield_geq_bicub(g,RPZ(1:2:N-2),RPZ(3:2:N),nowarn);
if ierr_Bas == 1
    if ~nowarn
        warning('AS bfield error in fl_derivs_phi_g3d')
    end
    ierr = 1; df = [];
    return;
end    

% NA part
if abs(max(rmp.current)) > 1e-8
   [Br,Bphi,Bz]=bfield_bs_cyl(RPZ(1:3:N-2),RPZ(2:3:N-1),RPZ(3:3:N),rmp.coil,rmp.current,nowarn);
   Bout.br = Bout.br + Br;
   Bout.bphi = Bout.bphi + Bphi;
   Bout.bz = Bout.bz + Bz;
end

Btot = sqrt(Bout.br.^2 + Bout.bz.^2 + Bout.bphi.^2);
df(1:3:N-2) = Bout.br./Btot; % dR/dl
df(2:3:N-1) = Bout.bphi./(RPZ(1:3:N-2).'.*Btot); % dphi/dl
df(3:3:N)   = Bout.bz./Btot; % dZ/dl
ierr = 0;



