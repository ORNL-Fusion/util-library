function g = correct_gfile_wb_to_wb_per_radian(g)
% Convert a gfile whose flux and flux derivatives use Wb to Wb/rad.
% Call this only after confirming the source convention. The Bpol check is
% advisory and does not prevent an explicitly requested conversion.

if isfield(g,'two_pi_corrected')
    error('This gfile has already been marked as two_pi_corrected')
end

% Check that the flux grid, header values, and derivative profiles can all
% be converted together. A partial conversion would leave the gfile inconsistent.
if ~isequal(size(g.psirz),[g.mw,g.mh]) || numel(g.ffprim) ~= g.mw || numel(g.pprime) ~= g.mw
    error('Gfile flux grid or derivative profiles have unexpected sizes')
end
if ~isscalar(g.ssimag) || ~isscalar(g.ssibry) || g.ssimag == g.ssibry
    error('Gfile axis and boundary flux values must be distinct scalars')
end
if any(~isfinite(g.psirz(:))) || any(~isfinite([g.ssimag,g.ssibry,g.ffprim(:).',g.pprime(:).']))
    error('Gfile flux values and derivative profiles must be finite')
end

% A raw gfile needs interpolation coefficients before the advisory Bpol check.
if ~isfield(g,'bicub_coeffs_inv')
    g = postprocess_gfile(g);
end
Diagnostic = diagnose_gfile_psi_units(g);
if ~Diagnostic.available
    fprintf(2,'WARNING: Could not check gfile flux units before conversion: %s. Converting anyway.\n',Diagnostic.reason);
elseif strcmp(Diagnostic.inferredUnits,'Wb/rad')
    fprintf(2,'WARNING: Bpol diagnostic favors Wb/rad over Wb; converting anyway as requested.\n');
end

% Flux values divide by 2*pi; derivatives with respect to flux multiply by
% 2*pi. Profiles tabulated on normalized flux do not change.
twoPi = 2*pi;
g.psirz = g.psirz/twoPi;
g.ssimag = g.ssimag/twoPi;
g.ssibry = g.ssibry/twoPi;
g.ffprim = g.ffprim*twoPi;
g.pprime = g.pprime*twoPi;
g.two_pi_corrected = true;

% Rebuild coefficients and diagnostics from the corrected flux grid.
g = postprocess_gfile(g);
end
