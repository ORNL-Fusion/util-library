function Diagnostic = diagnose_gfile_psi_units(g)
% Compare the two common gfile psi conventions with a rough Ip-based Bpol estimate.
% The circular-current estimate is only a diagnostic: shaping and external
% fields can change Bpol enough that its preferred units are not definitive.

Diagnostic = struct('available',false,'reason','','inferredUnits','','bpolEstimate',nan,'bpolWbPerRad',nan,'bpolWb',nan);

% A missing or zero plasma current cannot provide an independent field scale.
if isempty(g.cpasma) || ~isscalar(g.cpasma) || ~isfinite(g.cpasma) || abs(g.cpasma) < 1e-6
    Diagnostic.reason = 'Gfile does not contain a usable nonzero plasma current';
    return
end
if ~isscalar(g.rmaxis) || ~isscalar(g.zmaxis) || ~all(isfinite([g.rmaxis,g.zmaxis]))
    Diagnostic.reason = 'Gfile does not contain a usable magnetic axis';
    return
end
if ~isnumeric(g.bdry) || size(g.bdry,1) ~= 2 || size(g.bdry,2) < 2 || any(~isfinite(g.bdry(:)))
    Diagnostic.reason = 'Gfile does not contain a usable plasma boundary';
    return
end

% Find the outboard LCFS crossing at the magnetic-axis height. If the contour
% has no such crossing, use its largest-R point as in the SOLPS-routines logic.
[r,z] = outboard_boundary_at_axis(g.bdry,g.zmaxis);
rminor = hypot(r-g.rmaxis,z-g.zmaxis);
if r <= 0 || rminor <= 0
    Diagnostic.reason = 'Outboard boundary point is not usable';
    return
end

% Raw gfile psi is unscaled here. A Wb grid gives a field smaller by 2*pi
% than the same numerical gradient interpreted as Wb/rad.
try
    [dpsidr,dpsidz,ierr] = calc_psi_derivs(g,r,z,1);
catch ME
    Diagnostic.reason = ['Could not evaluate psi derivatives: ',ME.message];
    return
end
if ierr || ~isfinite(dpsidr) || ~isfinite(dpsidz)
    Diagnostic.reason = 'Could not evaluate psi derivatives at the outboard boundary';
    return
end
gradPsi = hypot(dpsidr,dpsidz);
bpolWbPerRad = gradPsi/r;
bpolWb = bpolWbPerRad/(2*pi);
mu0 = 4*pi*1e-7;
bpolEstimate = mu0*abs(g.cpasma)/(2*pi*rminor);
if ~all(isfinite([bpolWbPerRad,bpolWb,bpolEstimate])) || min([bpolWbPerRad,bpolWb,bpolEstimate]) <= 0
    Diagnostic.reason = 'Could not evaluate finite positive Bpol values';
    return
end

% Log-ratio distance treats overestimates and underestimates symmetrically.
scoreWbPerRad = abs(log(bpolWbPerRad/bpolEstimate));
scoreWb = abs(log(bpolWb/bpolEstimate));
Diagnostic.available = true;
Diagnostic.bpolEstimate = bpolEstimate;
Diagnostic.bpolWbPerRad = bpolWbPerRad;
Diagnostic.bpolWb = bpolWb;
if scoreWb < scoreWbPerRad
    Diagnostic.inferredUnits = 'Wb';
else
    Diagnostic.inferredUnits = 'Wb/rad';
end
end

function [r,z] = outboard_boundary_at_axis(bdry,zAxis)
% Return the largest-R intersection with zAxis, including the closing edge.
rCandidates = [];
nBoundary = size(bdry,2);
for iPoint = 1:nBoundary
    nextPoint = mod(iPoint,nBoundary)+1;
    r0 = bdry(1,iPoint);
    z0 = bdry(2,iPoint);
    r1 = bdry(1,nextPoint);
    z1 = bdry(2,nextPoint);
    if z0 == z1
        if z0 == zAxis
            rCandidates = [rCandidates,r0,r1];
        end
    elseif min(z0,z1) <= zAxis && zAxis <= max(z0,z1)
        fraction = (zAxis-z0)/(z1-z0);
        rCandidates(end+1) = r0+fraction*(r1-r0);
    end
end
if isempty(rCandidates)
    [r,iPoint] = max(bdry(1,:));
    z = bdry(2,iPoint);
else
    r = max(rCandidates);
    z = zAxis;
end
end
