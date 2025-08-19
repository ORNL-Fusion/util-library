function [df,ierr] = fl_derivs_dlp_gfile(RZ,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end
N = length(RZ);
Neq = 2;

% Axisym part
[Bout,ierr_Bas] = bfield_geq_bicub(bfield.g,RZ(1:Neq:N-1),RZ(2:Neq:N),nowarn);
if ierr_Bas == 1
    if ~nowarn
        warning('AS bfield error in fl_derivs_dlp_gfile')
    end
    ierr = 1; df = [];
    return;
end    

Bpol = max(hypot(Bout.br,Bout.bz),1e-12);
df = zeros(1,N);
df(1:Neq:N-1) = Bout.br./Bpol; % dR/dlp
df(2:Neq:N)   = Bout.bz./Bpol; % dZ/dlp
% dphi/dlp would be Bout.Bphi./(R.*Bpol);
ierr = 0;



