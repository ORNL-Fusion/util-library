function [df,ierr] = fl_derivs_dphi_just_coils(phi,RZ,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end
N = length(RZ);

if isfield(bfield,'bounds')
    switch bfield.bounds.type
        case 'box'
            if any(RZ(1:2:N-1) < bfield.bounds.Rlims(1) | RZ(1:2:N-1) > bfield.bounds.Rlims(2) | RZ(2:2:N) < bfield.bounds.Zlims(1) | RZ(2:2:N) > bfield.bounds.Zlims(2))
                if ~nowarn
                    fprintf('RZ out of defined bfield.bounds.[R,Z]lim\n')
                end
                df = [];
                ierr = 1;
                return
            end
        case 'ves'
            ves_cut = cut_W7X_vessel(bfield.bounds.ves,phi);
            if any(~inpolygon(RZ(1:2:N-1),RZ(2:2:N),ves_cut.r,ves_cut.z))
                if ~nowarn
                    fprintf('RZ out of defined bfield.bounds.[R,Z]lim\n')
                end
                df = [];
                ierr = 1;
                return
            end
        otherwise
            error('Did not recognize bfield.bounds.type')
    end
end

[Br,Bphi,Bz]=bfield_bs_cyl(RZ(1:2:N-1),phi*ones(size(RZ(1:2:N-1))),RZ(2:2:N),bfield.coil,bfield.current,nowarn);
df(1:2:N-1) = RZ(1:2:N-1).'.*Br./Bphi;
df(2:2:N)   = RZ(1:2:N-1).'.*Bz./Bphi;
ierr = 0;