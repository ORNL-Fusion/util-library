function [df,ierr] = fl_derivs_dphi_g3d(phi,RZ,g,rmp,nowarn)
if nargin < 5
    nowarn = 0;
end
N = length(RZ);

% Axisym part
[Bout,ierr_Bas] = bfield_geq_bicub(g,RZ(1:2:N-1),RZ(2:2:N),nowarn);
if ierr_Bas == 1
    if ~nowarn
        warning('AS bfield error in fl_derivs_phi_g3d')
    end
    ierr = 1; df = [];
    return;
end    

% NA part
if abs(max(rmp.current)) > 1e-8
   [Br,Bphi,Bz]=bfield_bs_cyl(RZ(1:2:N-1),phi,RZ(2:2:N),rmp.coil,rmp.current,nowarn);
   Bout.br = Bout.br + Br;
   Bout.bphi = Bout.bphi + Bphi;
   Bout.bz = Bout.bz + Bz;
end
   
df(1:2:N-1) = RZ(1:2:N-1).'.*Bout.br./Bout.bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bout.bz./Bout.bphi;
ierr = 0;