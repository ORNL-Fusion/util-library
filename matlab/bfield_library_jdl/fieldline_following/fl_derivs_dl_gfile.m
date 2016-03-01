function [df,ierr] = fl_derivs_dl_gfile(RPZ,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end
N = length(RPZ);

% Axisym part
[Bout,ierr_Bas] = bfield_geq_bicub(bfield.g,RPZ(1:2:N-2),RPZ(3:2:N),nowarn);
if ierr_Bas == 1
    if ~nowarn
        warning('AS bfield error in fl_derivs_dl_gfile')
    end
    ierr = 1; df = [];
    return;
end    

Btot = sqrt(Bout.br.^2 + Bout.bz.^2 + Bout.bphi.^2);
df(1:3:N-2) = Bout.br./Btot; % dR/dl
df(2:3:N-1) = Bout.bphi./(RPZ(1:3:N-2).'.*Btot); % dphi/dl
df(3:3:N)   = Bout.bz./Btot; % dZ/dl
ierr = 0;



