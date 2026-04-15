function [df,ierr] = fl_derivs_dphi_general(phi,RZ,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end

N = length(RZ);
[Bout,ierr] = bfield_general_rzphi(RZ(1:2:N-1),RZ(2:2:N),phi,bfield,nowarn);
if ierr == 1
    if ~nowarn
        warning('bfield error in fl_derivs_dphi_general')
    end
    df = [];
    return;
end

df = zeros(1,N);
df(1:2:N-1) = RZ(1:2:N-1).*Bout.br./Bout.bphi;
df(2:2:N)   = RZ(1:2:N-1).*Bout.bz./Bout.bphi;
