function Geo = lookup_flux_surface_geometry(g,R,Z)
% Interpolate cached closed-core geometry at paired R,Z points.
% Build the table with g = postprocess_gfile_geometry(g). R and Z must
% have matching sizes, or one may be scalar. Geo retains the expanded shape.
% Points must lie strictly inside the sampled psiN interval and on the
% closed-core branch. No extrapolation or endpoint
% clamping is used. Geometry fields are NaN when Geo.iErr = 1; errorId
% distinguishes InvalidPoint, OutsideGrid, OutsideTable, WrongBranch and
% InvalidInterpolation. Ambiguous branch checks can also return NoContour,
% NoMatchingContour or OpenContour. Stale or incomplete tables cause an error.
% Returns the OMP fields, both flux radii and their metrics. Contour arrays
% remain empty; compute_flux_surface_geometry returns query-level contours.

if ~isnumeric(R) || ~isnumeric(Z) || isempty(R) || isempty(Z) || ~isreal(R) || ~isreal(Z)
    error('lookup_flux_surface_geometry:InvalidPoint','R and Z must be real numeric arrays')
end
if isscalar(R) && ~isscalar(Z)
    R = repmat(R,size(Z));
elseif isscalar(Z) && ~isscalar(R)
    Z = repmat(Z,size(R));
elseif ~isequal(size(R),size(Z))
    error('lookup_flux_surface_geometry:InvalidPoint','R and Z must have matching sizes, or one must be scalar')
end

% Validate the cached equilibrium before using its bicubic flux interpolant.
% A NaN query supplies the same empty output template as the shared lookup.
Empty = lookup_flux_surface_table(g,NaN);
Table = g.fluxSurfaceGeometry;
rPoints = R(:).';
zPoints = Z(:).';
finitePoints = isfinite(rPoints) & isfinite(zPoints);
psiN = nan(size(rPoints));
outsideGrid = false(size(rPoints));
if any(finitePoints)
    [psiNFinite,~,ierrFinite] = calc_psiN(g,rPoints(finitePoints),zPoints(finitePoints),1);
    psiN(finitePoints) = psiNFinite;
    outsideGrid(finitePoints) = ierrFinite | ~isfinite(psiNFinite);
end
Geo = reshape(lookup_flux_surface_table(g,psiN),size(R));

% A flux value alone cannot distinguish the core from private-flux loops.
% Neighboring stored contours usually establish the closed-core branch.
% When polygon discretization exceeds the gap between nearby surfaces,
% resolve the branch from a query-level contour on the same refined grid.
BranchGrid = struct();
for iPoint = 1:numel(rPoints)
    errorId = '';
    if ~finitePoints(iPoint)
        errorId = 'InvalidPoint';
    elseif outsideGrid(iPoint)
        errorId = 'OutsideGrid';
    elseif Geo(iPoint).iErr == 0
        iInner = find(Table.psiN < psiN(iPoint),1,'last');
        iOuter = find(Table.psiN > psiN(iPoint),1,'first');
        insideOuter = inpolygon(rPoints(iPoint),zPoints(iPoint),Table.contourR{iOuter},Table.contourZ{iOuter});
        insideInner = inpolygon(rPoints(iPoint),zPoints(iPoint),Table.contourR{iInner},Table.contourZ{iInner});
        if ~insideOuter || insideInner
            if isempty(fieldnames(BranchGrid))
                [BranchGrid.r,BranchGrid.z,psiGrid] = refine_psiN(g.r,g.z,Table.refinement,g);
                BranchGrid.psiN = psiGrid.';
                BranchGrid.spacing = hypot(BranchGrid.r(2)-BranchGrid.r(1),BranchGrid.z(2)-BranchGrid.z(1));
            end
            errorId = query_core_branch(g,rPoints(iPoint),zPoints(iPoint),psiN(iPoint),BranchGrid);
        end
    end
    if ~isempty(errorId)
        Geo(iPoint) = Empty;
        Geo(iPoint).psiN = psiN(iPoint);
        Geo(iPoint).errorId = errorId;
    end
end
end

function errorId = query_core_branch(g,R,Z,psiN,Grid)
% Choose the closest contour component, require closure, and check that
% it encloses the axis. The matching and closure tolerances use the same
% grid-spacing convention as compute_flux_surface_geometry; area integrals
% are unnecessary for branch selection. Private-flux loops remain rejected.
errorId = 'NoContour';
contours = contourc(Grid.r,Grid.z,Grid.psiN,[psiN,psiN]);
bestDistance = inf;
bestContour = [];
iColumn = 1;
while iColumn <= size(contours,2)
    nVertices = contours(2,iColumn);
    points = contours(:,iColumn+1:iColumn+nVertices);
    iColumn = iColumn+nVertices+1;
    if nVertices < 2 || any(~isfinite(points(:)))
        continue
    end
    startPoints = points(:,1:end-1);
    segments = points(:,2:end)-startPoints;
    lengthSquared = sum(segments.^2,1);
    lengthSquared(lengthSquared == 0) = 1;
    fraction = ((R-startPoints(1,:)).*segments(1,:)+(Z-startPoints(2,:)).*segments(2,:))./lengthSquared;
    closest = startPoints+segments.*max(0,min(1,fraction));
    distance = min(hypot(R-closest(1,:),Z-closest(2,:)));
    if distance < bestDistance
        bestDistance = distance;
        bestContour = points;
    end
end
if isempty(bestContour)
    return
end
if bestDistance > 0.25*Grid.spacing
    errorId = 'NoMatchingContour';
elseif size(bestContour,2) < 4 || hypot(bestContour(1,1)-bestContour(1,end),bestContour(2,1)-bestContour(2,end)) > 1e-6*Grid.spacing
    errorId = 'OpenContour';
elseif ~inpolygon(g.rmaxis,g.zmaxis,bestContour(1,:),bestContour(2,:))
    errorId = 'WrongBranch';
else
    errorId = '';
end
end
