function [df,ierr] = fl_derivs_dphi_gfile(RZ,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end
N = length(RZ);

% Axisym part
Neq = 2;
[Bout,ierr] = bfield_geq_bicub(bfield.g,RZ(1:Neq:N-1),RZ(2:Neq:N),nowarn);
if ierr == 1
    if ~nowarn
        warning('AS bfield error in fl_derivs_dphi_gfile')
    end
    ierr = 1; df = [];
    return;
end    

df = zeros(1,N);
df(1:2:N-1) = RZ(1:2:N-1).*Bout.br./Bout.bphi;
df(2:2:N)   = RZ(1:2:N-1).*Bout.bz./Bout.bphi;
ierr = 0;