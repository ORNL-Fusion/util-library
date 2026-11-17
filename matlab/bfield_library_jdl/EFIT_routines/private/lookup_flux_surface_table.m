function Geo = lookup_flux_surface_table(g,psiN)
% Shared interpolation for geometry and radial-label queries on core surfaces.
% Coordinate values and their contour-derived slopes define cubic Hermite
% interpolants in psiN. Differentiating those same interpolants keeps label
% derivatives and metric transformations consistent between table samples.
if ~isnumeric(psiN) || ~isreal(psiN) || isempty(psiN)
    error('lookup_flux_surface_geometry:InvalidPoint','psiN must be a nonempty real numeric array')
end
Table = validated_table(g);
Geo = repmat(empty_geometry(Table),size(psiN));
p = psiN(:).';
finitePoints = isfinite(p);
valid = finitePoints & p > Table.psiN(1) & p < Table.psiN(end);
for iPoint = 1:numel(p)
    Geo(iPoint).psiN = p(iPoint);
    if ~finitePoints(iPoint)
        Geo(iPoint).errorId = 'InvalidPoint';
    elseif ~valid(iPoint)
        Geo(iPoint).errorId = 'OutsideTable';
    end
end
iValid = find(valid);
if isempty(iValid)
    return
end

% The base coordinate is p=psiN. Vp=dV/dp and Kp=Vp*<|grad p|^2>
% are independent of which radial label is requested. Interpolate these
% once, then transform with the derivatives of the coordinate interpolants.
vp = interp1(Table.psiN,Table.dVdx./Table.dpsiNdx,p(iValid),'pchip');
kp = interp1(Table.psiN,Table.transportMetric.*Table.dpsiNdx,p(iValid),'pchip');
areaRZ = interp1(Table.psiN,Table.areaRZ,p(iValid),'pchip');
volume = interp1(Table.psiN,Table.volume,p(iValid),'pchip');
centroidR = interp1(Table.psiN,Table.centroidR,p(iValid),'pchip');
dxdpSamples = 1./Table.dpsiNdx;
dRadiusSquaredDp = 2*Table.rSaarelma.*Table.drSaarelmaDx./Table.dpsiNdx;
for iQuery = 1:numel(iValid)
    iPoint = iValid(iQuery);
    [x,dxdp] = coordinate_value_and_slope(Table.psiN,Table.x,dxdpSamples,p(iPoint));
    [phiTor,dPhiTorDpsiN] = coordinate_value_and_slope(Table.psiN,Table.phiTor,Table.dPhiTorDpsiN,p(iPoint));
    [radiusSquared,dRadiusSquared] = coordinate_value_and_slope(Table.psiN,Table.rSaarelma.^2,dRadiusSquaredDp,p(iPoint));
    positiveValues = [dxdp,phiTor,dPhiTorDpsiN,radiusSquared,dRadiusSquared,vp(iQuery),kp(iQuery),areaRZ(iQuery),volume(iQuery),centroidR(iQuery)];
    if ~isfinite(x) || any(~isfinite(positiveValues) | positiveValues <= 0)
        Geo(iPoint).errorId = 'InvalidInterpolation';
        continue
    end

    % Interpolate squared radii through their defining enclosed integrals.
    % This preserves rTor^2=phiTor/(pi*B0) and avoids an independent radius
    % fit disagreeing with its flux or its derivative.
    rTor = sqrt(phiTor/(pi*Table.b0Reference));
    rSaarelma = sqrt(radiusSquared);
    drTorDp = dPhiTorDpsiN/(2*pi*Table.b0Reference*rTor);
    drSaarelmaDp = dRadiusSquared/(2*rSaarelma);
    drTorDx = drTorDp/dxdp;
    drSaarelmaDx = drSaarelmaDp/dxdp;
    dVdx = vp(iQuery)/dxdp;
    transportMetric = kp(iQuery)*dxdp;
    gradXSqAvg = transportMetric/dVdx;
    dVdrTor = vp(iQuery)/drTorDp;
    dVdrSaarelma = vp(iQuery)/drSaarelmaDp;
    transportMetricRTor = kp(iQuery)*drTorDp;
    transportMetricRSaarelma = kp(iQuery)*drSaarelmaDp;
    gradRTorSqAvg = transportMetricRTor/dVdrTor;
    gradRSaarelmaSqAvg = transportMetricRSaarelma/dVdrSaarelma;
    transformedValues = [rTor,rSaarelma,1/dxdp,drTorDx,drSaarelmaDx,dVdx,transportMetric,gradXSqAvg,dVdrTor,dVdrSaarelma,transportMetricRTor,transportMetricRSaarelma,gradRTorSqAvg,gradRSaarelmaSqAvg];
    if any(~isfinite(transformedValues) | transformedValues <= 0)
        Geo(iPoint).errorId = 'InvalidInterpolation';
        continue
    end

    Geo(iPoint).x = x;
    Geo(iPoint).areaRZ = areaRZ(iQuery);
    Geo(iPoint).volume = volume(iQuery);
    Geo(iPoint).centroidR = centroidR(iQuery);
    Geo(iPoint).dpsiNdx = 1/dxdp;
    Geo(iPoint).dVdx = dVdx;
    Geo(iPoint).gradXSqAvg = gradXSqAvg;
    Geo(iPoint).transportMetric = transportMetric;
    Geo(iPoint).phiTor = phiTor;
    Geo(iPoint).rTor = rTor;
    Geo(iPoint).rSaarelma = rSaarelma;
    Geo(iPoint).dPhiTorDpsiN = dPhiTorDpsiN;
    Geo(iPoint).drTorDx = drTorDx;
    Geo(iPoint).drSaarelmaDx = drSaarelmaDx;
    Geo(iPoint).dVdrTor = dVdrTor;
    Geo(iPoint).dVdrSaarelma = dVdrSaarelma;
    Geo(iPoint).gradRTorSqAvg = gradRTorSqAvg;
    Geo(iPoint).gradRSaarelmaSqAvg = gradRSaarelmaSqAvg;
    Geo(iPoint).transportMetricRTor = transportMetricRTor;
    Geo(iPoint).transportMetricRSaarelma = transportMetricRSaarelma;
    Geo(iPoint).iErr = 0;
    Geo(iPoint).errorId = '';
end
end

function [value,slope] = coordinate_value_and_slope(pNodes,values,slopes,p)
% The cubic matches both physical values and first derivatives at each end.
% Strictly interior queries always have a bracketing interval. Invalid
% interpolated slopes are reported by the caller rather than clipped.
iLeft = find(pNodes < p,1,'last');
h = pNodes(iLeft+1)-pNodes(iLeft);
delta = (values(iLeft+1)-values(iLeft))/h;
c2 = (3*delta-2*slopes(iLeft)-slopes(iLeft+1))/h;
c3 = (slopes(iLeft)+slopes(iLeft+1)-2*delta)/h^2;
offset = p-pNodes(iLeft);
value = values(iLeft)+offset*(slopes(iLeft)+offset*(c2+offset*c3));
slope = slopes(iLeft)+offset*(2*c2+3*offset*c3);
end

function Table = validated_table(g)
% All entry points require the same versioned geometry and source snapshot.
% Comparing raw equilibrium inputs catches field or geometry edits made
% after the table was built, without modifying g or writing a cache.
if ~isstruct(g) || ~isscalar(g) || ~isfield(g,'fluxSurfaceGeometry')
    error('lookup_flux_surface_geometry:NoTable','Run postprocess_gfile_geometry before lookup')
end
Table = g.fluxSurfaceGeometry;
sampleFields = {'x','areaRZ','volume','dVdx','gradXSqAvg','dpsiNdx','transportMetric','centroidR','phiTor','rTor','rSaarelma','dPhiTorDpsiN','drTorDx','drSaarelmaDx','dVdrTor','dVdrSaarelma','gradRTorSqAvg','gradRSaarelmaSqAvg','transportMetricRTor','transportMetricRSaarelma'};
metadataFields = {'version','Source','psiN','contourR','contourZ','refinement','zOMP','xDefinition','b0Reference','r0Reference','toroidalFieldSign','phiTorDefinition','rTorDefinition','rSaarelmaDefinition'};
if ~isstruct(Table) || ~isscalar(Table) || ~all(isfield(Table,[sampleFields,metadataFields])) || ~isequal(Table.version,2)
    error('lookup_flux_surface_geometry:InvalidTable','Rebuild the geometry table with postprocess_gfile_geometry')
end
if ~isnumeric(Table.psiN) || ~isreal(Table.psiN) || ~isvector(Table.psiN) || numel(Table.psiN) < 2 || any(~isfinite(Table.psiN(:))) || any(diff(Table.psiN(:)) <= 0) || any(Table.psiN(:) <= 0 | Table.psiN(:) >= 1)
    error('lookup_flux_surface_geometry:InvalidTable','The table must contain ordered interior flux surfaces')
end
Table.psiN = Table.psiN(:).';
nSamples = numel(Table.psiN);
for iField = 1:numel(sampleFields)
    fieldName = sampleFields{iField};
    values = Table.(fieldName);
    if ~isnumeric(values) || ~isreal(values) || ~isvector(values) || numel(values) ~= nSamples || any(~isfinite(values(:)))
        error('lookup_flux_surface_geometry:InvalidTable','Every geometry field must have one finite value per sampled surface')
    end
    Table.(fieldName) = values(:).';
end
if any(Table.dVdx <= 0 | Table.dpsiNdx <= 0 | Table.transportMetric <= 0 | Table.phiTor <= 0 | Table.rTor <= 0 | Table.rSaarelma <= 0 | Table.dPhiTorDpsiN <= 0 | Table.drTorDx <= 0 | Table.drSaarelmaDx <= 0) || any(diff(Table.x) <= 0) || any(diff(Table.phiTor) <= 0) || any(diff(Table.rSaarelma) <= 0)
    error('lookup_flux_surface_geometry:InvalidTable','The cached radial mappings must be positive and monotonic')
end
if ~iscell(Table.contourR) || ~iscell(Table.contourZ) || numel(Table.contourR) ~= nSamples || numel(Table.contourZ) ~= nSamples || ~isscalar(Table.b0Reference) || ~isfinite(Table.b0Reference) || Table.b0Reference <= 0 || ~isscalar(Table.r0Reference) || ~isfinite(Table.r0Reference) || Table.r0Reference <= 0
    error('lookup_flux_surface_geometry:InvalidTable','The cached contours or reference field are incomplete')
end
sourceFields = {'psirz','ssimag','ssibry','rmaxis','zmaxis','r','z','ip_sign','pn','fpol','rzero','bcentr'};
if ~isstruct(Table.Source) || ~isscalar(Table.Source) || ~all(isfield(Table.Source,sourceFields)) || ~all(isfield(g,sourceFields))
    error('lookup_flux_surface_geometry:InvalidTable','The geometry source snapshot is incomplete')
end
for iField = 1:numel(sourceFields)
    fieldName = sourceFields{iField};
    if ~isequaln(Table.Source.(fieldName),g.(fieldName))
        error('lookup_flux_surface_geometry:StaleTable','Equilibrium inputs changed; run postprocess_gfile_geometry again')
    end
end
end

function Geo = empty_geometry(Table)
% Physical quantities start unavailable; reference metadata and definitions
% remain available even when an individual query is outside the table.
Geo = struct('psiN',NaN,'x',NaN,'areaRZ',NaN,'volume',NaN,'dVdx',NaN,'gradXSqAvg',NaN,'dpsiNdx',NaN,'transportMetric',NaN,'centroidR',NaN,'contourR',[],'contourZ',[],'refinement',Table.refinement,'zOMP',Table.zOMP,'xDefinition',Table.xDefinition,'iErr',1,'errorId','');
Geo.phiTor = NaN;
Geo.rTor = NaN;
Geo.rSaarelma = NaN;
Geo.dPhiTorDpsiN = NaN;
Geo.drTorDx = NaN;
Geo.drSaarelmaDx = NaN;
Geo.dVdrTor = NaN;
Geo.dVdrSaarelma = NaN;
Geo.gradRTorSqAvg = NaN;
Geo.gradRSaarelmaSqAvg = NaN;
Geo.transportMetricRTor = NaN;
Geo.transportMetricRSaarelma = NaN;
Geo.b0Reference = Table.b0Reference;
Geo.r0Reference = Table.r0Reference;
Geo.toroidalFieldSign = Table.toroidalFieldSign;
Geo.phiTorDefinition = Table.phiTorDefinition;
Geo.rTorDefinition = Table.rTorDefinition;
Geo.rSaarelmaDefinition = Table.rSaarelmaDefinition;
end
