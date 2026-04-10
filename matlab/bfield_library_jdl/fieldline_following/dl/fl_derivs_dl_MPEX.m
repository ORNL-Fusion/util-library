function [df,ierr] = fl_derivs_dl_MPEX(RPZ,bfield,nowarn)
if nargin < 3
    nowarn = 0;
end
N = length(RPZ);

r = RPZ(1:3:N-2);
z = RPZ(3:3:N);

if isfield(bfield,'dcoil_cutoff') && ~isempty(bfield.dcoil_cutoff) && bfield.dcoil_cutoff > 0
    rwind = reshape(bfield.coil.rwind,1,[]);
    dwind = reshape(bfield.coil.zwind,1,[]);
    dcoil = sqrt((r(:) - rwind).^2 + (z(:) - dwind).^2);
    if any(dcoil(:) < bfield.dcoil_cutoff)
        if ~nowarn
            warning('Field-line integration stopped near a coil filament.')
        end
        df = NaN(size(RPZ));
        ierr = 1;
        return
    end
end

[Br,Bz] = bfield_circular_coils_vectorized_zconst(bfield.coil,bfield.current,r,z);
Bphi = zeros(size(Br));
Btot = sqrt(Br.^2 + Bz.^2);
if any(~isfinite(Btot)) || any(Btot <= 0)
    if ~nowarn
        warning('Field-line integration encountered invalid total field magnitude.')
    end
    df = NaN(size(RPZ));
    ierr = 1;
    return
end
	
df(1:3:N-2) = Br./Btot; % dR/dl
df(2:3:N-1) = Bphi./(r.'.*Btot); % dphi/dl
df(3:3:N)   = Bz./Btot; % dZ/dl
ierr = 0;
