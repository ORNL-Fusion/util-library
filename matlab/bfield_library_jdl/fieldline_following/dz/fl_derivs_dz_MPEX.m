function [df,ierr] = fl_derivs_dz_MPEX(Z,RP,bfield,nowarn)
if nargin < 4
    nowarn = 0;
end

N = length(RP);
r = RP(1:2:N-1);

if isfield(bfield,'dcoil_cutoff') && ~isempty(bfield.dcoil_cutoff) && bfield.dcoil_cutoff > 0
    rwind = reshape(bfield.coil.rwind,1,[]);
    dwind = reshape(bfield.coil.zwind,1,[]);
    dcoil = sqrt((r(:) - rwind).^2 + (Z - dwind).^2);
    if any(dcoil(:) < bfield.dcoil_cutoff)
        if ~nowarn
            warning('Field-line integration stopped near a coil filament.')
        end
        df = NaN(size(RP));
        ierr = 1;
        return
    end
end

[Br,Bz] = bfield_circular_coils_vectorized_zconst(bfield.coil,bfield.current,r,Z);
df(1:2:N-1) = Br./Bz; % dR/dz
df(2:2:N)   = zeros(size(Br)); % dphi/dz
ierr = 0;
