function Geo = lookup_flux_surface_geometry(g,R,Z)
% Interpolate a precomputed core geometry table at paired R,Z points.
% Build the table first with g = postprocess_gfile_geometry(g). Points must
% lie strictly within its sampled psiN range and between the neighboring
% closed core contours. Exact query-level contours are not returned; use
% compute_flux_surface_geometry for those. Failures set Geo.iErr = 1.

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
if ~isstruct(g) || ~isfield(g,'fluxSurfaceGeometry')
    error('lookup_flux_surface_geometry:NoTable','Run postprocess_gfile_geometry before lookup')
end
Table = g.fluxSurfaceGeometry;
requiredFields = {'psiN','x','areaRZ','volume','dVdx','gradXSqAvg','dpsiNdx','transportMetric','centroidR','contourR','contourZ','refinement','zOMP','xDefinition'};
if ~isstruct(Table) || ~all(isfield(Table,requiredFields)) || numel(Table.psiN) < 2 || any(diff(Table.psiN) <= 0)
    error('lookup_flux_surface_geometry:InvalidTable','The gfile geometry table is incomplete or unordered')
end

% Evaluate normalized flux only for finite coordinates. Query-level
% contours are replaced by containment between neighboring core loops.
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

Geo = repmat(empty_geometry(Table),size(R));
validLookup = false(size(rPoints));
for iPoint = 1:numel(rPoints)
    Geo(iPoint).psiN = psiN(iPoint);
    if ~finitePoints(iPoint)
        Geo(iPoint).errorId = 'InvalidPoint';
    elseif outsideGrid(iPoint)
        Geo(iPoint).errorId = 'OutsideGrid';
    elseif psiN(iPoint) <= Table.psiN(1) || psiN(iPoint) >= Table.psiN(end)
        Geo(iPoint).errorId = 'OutsideTable';
    else
        iInner = find(Table.psiN < psiN(iPoint),1,'last');
        iOuter = find(Table.psiN > psiN(iPoint),1,'first');
        insideOuter = inpolygon(rPoints(iPoint),zPoints(iPoint),Table.contourR{iOuter},Table.contourZ{iOuter});
        insideInner = inpolygon(rPoints(iPoint),zPoints(iPoint),Table.contourR{iInner},Table.contourZ{iInner});
        if ~insideOuter || insideInner
            Geo(iPoint).errorId = 'WrongBranch';
        else
            validLookup(iPoint) = true;
            Geo(iPoint).iErr = 0;
        end
    end
end

% Interpolate all accepted points together. Interpolating stored geometry
% avoids extracting a new contour or evaluating a flux-surface integral.
iValid = find(validLookup);
if isempty(iValid)
    return
end
fields = {'x','areaRZ','volume','dVdx','dpsiNdx','transportMetric','centroidR'};
for iField = 1:numel(fields)
    fieldName = fields{iField};
    values = interp1(Table.psiN,Table.(fieldName),psiN(iValid),'pchip');
    for i = 1:numel(iValid)
        Geo(iValid(i)).(fieldName) = values(i);
    end
end
for i = 1:numel(iValid)
    iPoint = iValid(i);
    Geo(iPoint).gradXSqAvg = Geo(iPoint).transportMetric/Geo(iPoint).dVdx;
end
end

function Geo = empty_geometry(Table)
% Match the direct routine's scalar geometry fields and status convention.
Geo = struct('psiN',NaN,'x',NaN,'areaRZ',NaN,'volume',NaN,'dVdx',NaN,'gradXSqAvg',NaN,'dpsiNdx',NaN,'transportMetric',NaN,'centroidR',NaN,'contourR',[],'contourZ',[],'refinement',Table.refinement,'zOMP',Table.zOMP,'xDefinition',Table.xDefinition,'iErr',1,'errorId','');
end
