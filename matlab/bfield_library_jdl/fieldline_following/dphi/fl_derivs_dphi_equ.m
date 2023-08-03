function [df,ierr] = fl_derivs_dphi_equ(RZ,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end
N = length(RZ);

% Axisym part
[Bout,ierr] = bfield_equ_bicub(bfield.Equ,RZ(1:2:N-1),RZ(2:2:N),nowarn);
if ierr == 1
    if ~nowarn
        warning('AS bfield error in fl_derivs_dphi_equ')
    end
    ierr = 1; df = [];
    return;
end    

df(1:2:N-1) = RZ(1:2:N-1).*Bout.br./Bout.bphi;
df(2:2:N)   = RZ(1:2:N-1).*Bout.bz./Bout.bphi;
ierr = 0;