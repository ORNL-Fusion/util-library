function [df,ierr] = fl_derivs_dl_general(RPZ,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end

N = length(RPZ);
[Bout,ierr] = bfield_general_rzphi(RPZ(1:3:N-2),RPZ(3:3:N),RPZ(2:3:N-1),bfield,nowarn);
if ierr == 1
    if ~nowarn
        warning('bfield error in fl_derivs_dl_general')
    end
    df = [];
    return;
end

Btot = sqrt(Bout.br.^2 + Bout.bz.^2 + Bout.bphi.^2);
df = zeros(1,N);
df(1:3:N-2) = Bout.br./Btot;
df(2:3:N-1) = Bout.bphi./(RPZ(1:3:N-2).*Btot);
df(3:3:N)   = Bout.bz./Btot;
