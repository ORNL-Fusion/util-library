function Endpoints = compute_flux_surface_endpoints(g,refinement)
% Compute finite axis and separatrix flux labels without endpoint metrics.
% Usage: Endpoints = compute_flux_surface_endpoints(g,refinement)
% refinement is optional and follows compute_flux_surface_geometry.
% The axis flux and radii are zero. Separatrix values are limits of whole
% core area integrals at psiN=1-delta, delta=1e-3/2.^(0:5).
% For a regular X point I(delta)=L+a*delta*log(delta)+b*delta+... .
% I(delta)-4*I(delta/2)+4*I(delta/4) cancels both leading terms and also
% applies to a smooth boundary (a=0). Adjacent extrapolation windows and
% successive doubled contour grids must agree within 1e-4 relative for
% phiTor and rSaarelma^2. The reported change is a convergence indicator,
% not an error bound on the supplied equilibrium. Degenerate X points or
% unresolved topology can fail this check. Failed separatrix values are
% NaN; check separatrixIErr and separatrixErrorId. No cache is written.

if nargin < 2
    refinement = [];
end
Axis = compute_flux_surface_geometry(g,g.rmaxis,g.zmaxis,refinement);
Endpoints = struct('phiTorAxis',0,'rTorAxis',0,'rSaarelmaAxis',0,'phiTorSep',NaN,'rTorSep',NaN,'rSaarelmaSep',NaN,'separatrixRelativeChange',NaN,'separatrixRefinement',NaN,'separatrixIErr',1,'separatrixErrorId','SeparatrixNotConverged');
Endpoints.b0Reference = Axis.b0Reference;
Endpoints.r0Reference = Axis.r0Reference;
Endpoints.toroidalFieldSign = Axis.toroidalFieldSign;
factor = Axis.refinement;
levels = 1-1e-3./2.^(0:5);

% The outboard line selects the core surface at each level. All integrals
% use the resulting closed contour and include the entire axis region.
rLine = linspace(g.rmaxis,g.r(end),max(512,2*numel(g.r)));
psiLine = calc_psiN(g,rLine,g.zmaxis*ones(size(rLine)),1);
radii = nan(size(levels));
for iLevel = 1:numel(levels)
    iCross = find(psiLine(1:end-1) <= levels(iLevel) & psiLine(2:end) >= levels(iLevel),1,'first');
    if ~isempty(iCross)
        radii(iLevel) = fzero(@(rValue) calc_psiN(g,rValue,g.zmaxis,1)-levels(iLevel),rLine(iCross:iCross+1));
    end
end
previous = [];
previousChange = NaN;
for iGrid = 1:3
    Geo = compute_flux_surface_geometry(g,radii,g.zmaxis,factor);
    integrals = [[Geo.phiTor].',[Geo.rSaarelma].'.^2];
    bestChange = inf;
    estimate = [];
    for iWindow = 1:numel(levels)-3
        indices = iWindow:iWindow+3;
        values = integrals(indices,:);
        if any([Geo(indices).iErr]) || any(~isfinite(values(:)))
            continue
        end
        left = values(1,:)-4*values(2,:)+4*values(3,:);
        right = values(2,:)-4*values(3,:)+4*values(4,:);
        increments = diff(values,1,1);
        if any(right <= values(end,:)) || any(increments(:) <= 0)
            continue
        end
        change = max(abs(right-left)./right);
        if change < bestChange
            bestChange = change;
            estimate = right;
        end
    end
    if ~isempty(estimate)
        if ~isempty(previous)
            change = max([bestChange,previousChange,abs(estimate-previous)./estimate]);
            Endpoints.separatrixRelativeChange = change;
            if change <= 1e-4
                Endpoints.phiTorSep = estimate(1);
                Endpoints.rTorSep = sqrt(estimate(1)/(pi*Axis.b0Reference));
                Endpoints.rSaarelmaSep = sqrt(estimate(2));
                Endpoints.separatrixRefinement = factor;
                Endpoints.separatrixIErr = 0;
                Endpoints.separatrixErrorId = '';
                return
            end
        end
        previous = estimate;
        previousChange = bestChange;
    else
        previous = [];
    end
    Endpoints.separatrixRefinement = factor;
    factor = 2*factor;
end
end
