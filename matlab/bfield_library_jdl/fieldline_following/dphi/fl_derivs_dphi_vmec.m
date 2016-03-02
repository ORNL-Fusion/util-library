function [df,ierr] = fl_derivs_dphi_vmec(phi,RZ,bfield,nowarn)
if nargin < 5
    nowarn = 0;
end
N = length(RZ);

[Br,Bz,Bphi] = bcyl_vmec(RZ(1:2:N-1),(2:2:N),phi,bfield.wout);
ierr_Bvmec = 0;  % need to implement this
if ierr_Bvmec == 1
    if ~nowarn
        warning('bfield error in fl_derivs_dphi_vmec')
    end
    ierr = 1; df = [];
    return;
end
    
df(1:2:N-1) = RZ(1:2:N-1).'.*Br./Bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bz./Bphi;
ierr = 0;