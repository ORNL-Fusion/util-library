function [df,ierr] = fl_derivs_dphi_bspline(phi,RZ,g,bspline_data,nowarn)
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
if bspline_data.skip == 1
    Br = 0;
    Bz = 0;
    Bphi = 0;
else
    [Br,Bphi,Bz,ierr_B] = bfield_bs_Bspline(RZ(1:2:N-1),phi,RZ(2:2:N),bspline_data,nowarn);
%     ierr_B = 0;  % need to implement this
    if ierr_B == 1
        if ~nowarn
            warning('bfield error in fl_derivs_phi_bspline')
        end
        ierr = 1; df = [];
        return;
    end
end

Br = Bout.br + Br;
Bphi = Bout.bphi + Bphi;
Bz = Bout.bz + Bz;
    
df(1:2:N-1) = RZ(1:2:N-1).'.*Br./Bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bz./Bphi;
ierr = 0;