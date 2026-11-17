function Geo = compute_flux_surface_geometry(g,R,Z,refinement,zOMP)
% Compute geometry of core flux surfaces through (R,Z).
% R and Z must have matching sizes, or one may be scalar. Geo has the same
% size as the expanded inputs, with one contour and geometry per element.
% g must be a postprocessed gfile. The separatrix needs separate X-point
% handling. The optional refinement overrides automatic contour sampling.
% By default, use at least 256 points in one direction and spacing below
% 5 cm. Geo.refinement records the factor used. Expected point failures
% set Geo.iErr = 1 and Geo.errorId; geometry fields are NaN and contours empty.
% Geo.x is R_OMP-R_sep,OMP at zOMP (default g.zmaxis). Geo.dVdx and
% Geo.gradXSqAvg are dV/dx and <|grad x|^2> for this flux-surface label.
% Geo.transportMetric is their product for conservative transport terms.
% Geo.dpsiNdx is the OMP derivative. At the axis, gradXSqAvg is NaN;
% near the separatrix, check convergence with a larger refinement.
% phiTor is the positive magnitude of enclosed toroidal flux (Wb).
% rTor = sqrt(phiTor/(pi*b0Reference)); rSaarelma^2 = r0Reference/pi
% times the enclosed poloidal-area integral of 1/R. Both radii are in m.
% b0Reference = abs(g.bcentr*g.rzero/g.rmaxis), r0Reference = g.rmaxis.
% drTorDx and drSaarelmaDx transform the OMP metrics into dVdrTor,
% dVdrSaarelma, gradRTorSqAvg, gradRSaarelmaSqAvg and the corresponding
% transportMetricRTor and transportMetricRSaarelma products.
% At the axis the radii, flux, volume derivatives and transport products
% are zero; radius derivatives and averaged squared gradients are NaN.
% Area quadrature uses base-MATLAB polyshape triangulation (R2017b or later).

if ~isnumeric(R) || ~isnumeric(Z) || isempty(R) || isempty(Z) || ~isreal(R) || ~isreal(Z)
    error('compute_flux_surface_geometry:InvalidPoint','R and Z must be real numeric arrays')
end
if isscalar(R) && ~isscalar(Z)
    R = repmat(R,size(Z));
elseif isscalar(Z) && ~isscalar(R)
    Z = repmat(Z,size(R));
elseif ~isequal(size(R),size(Z))
    error('compute_flux_surface_geometry:InvalidPoint','R and Z must have matching sizes, or one must be scalar')
end

requiredFields = {'r','z','psirz','rmaxis','zmaxis','ssimag','ssibry','ip_sign','bicub_coeffs_inv','pn','fpol','rzero','bcentr'};
if ~isstruct(g) || ~all(isfield(g,requiredFields))
    error('compute_flux_surface_geometry:InvalidGfile','g must be a postprocessed gfile')
end
if g.ssimag == g.ssibry
    error('compute_flux_surface_geometry:InvalidGfile','Axis and boundary flux must differ')
end
% The reference field is the vacuum toroidal field at the axis major
% radius. A single sign of F is required for a monotonic toroidal-flux radius.
if ~isscalar(g.rmaxis) || ~isfinite(g.rmaxis) || g.rmaxis <= 0 || ~isscalar(g.rzero) || ~isfinite(g.rzero) || g.rzero <= 0 || ~isscalar(g.bcentr) || ~isfinite(g.bcentr) || g.bcentr == 0
    error('compute_flux_surface_geometry:InvalidReferenceField','rmaxis and rzero must be positive and bcentr must be finite and nonzero')
end
if ~isnumeric(g.pn) || ~isreal(g.pn) || ~isvector(g.pn) || numel(g.pn) < 2 || any(~isfinite(g.pn(:))) || any(diff(g.pn(:)) <= 0) || g.pn(1) ~= 0 || g.pn(end) ~= 1 || ~isnumeric(g.fpol) || ~isreal(g.fpol) || ~isvector(g.fpol) || numel(g.fpol) ~= numel(g.pn) || any(~isfinite(g.fpol(:))) || ~(all(g.fpol(:) > 0) || all(g.fpol(:) < 0))
    error('compute_flux_surface_geometry:InvalidToroidalField','pn must span 0 to 1 and fpol must be finite, nonzero and have a single sign')
end
Reference = struct('r0',g.rmaxis,'b0',abs(g.bcentr*g.rzero/g.rmaxis),'fieldSign',sign(g.fpol(1)));
if ~isfinite(Reference.b0) || Reference.b0 <= 0
    error('compute_flux_surface_geometry:InvalidReferenceField','The vacuum toroidal field at the axis must be finite and nonzero')
end
Fpol = pchip(g.pn,abs(g.fpol));

if nargin < 5
    zOMP = g.zmaxis;
end
if ~isscalar(zOMP) || ~isnumeric(zOMP) || ~isreal(zOMP) || ~isfinite(zOMP) || zOMP < g.z(1) || zOMP > g.z(end)
    error('compute_flux_surface_geometry:InvalidOmpHeight','zOMP must be a finite height inside the gfile grid')
end

% Choose the smallest integer factor satisfying the contour grid size and
% real-space spacing targets. refine_psiN uses refinement*nr by refinement*nz points.
if nargin < 4 || isempty(refinement)
    nr = numel(g.r);
    nz = numel(g.z);
    rExtent = abs(g.r(end)-g.r(1));
    zExtent = abs(g.z(end)-g.z(1));
    refinement = max(1,ceil(256/max(nr,nz)));
    while max([rExtent/(refinement*nr-1),zExtent/(refinement*nz-1)]) >= 0.05
        refinement = refinement+1;
    end
end
if ~isscalar(refinement) || ~isnumeric(refinement) || ~isreal(refinement) || ~isfinite(refinement) || refinement < 1 || refinement ~= floor(refinement)
    error('compute_flux_surface_geometry:InvalidRefinement','Refinement must be a positive integer')
end

% calc_psiN returns row vectors. Keep linear point order so Geo retains the
% shape of R and Z. Skip nonfinite coordinates before interpolation.
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
axisPoints = finitePoints & (rPoints == g.rmaxis) & (zPoints == g.zmaxis);
interiorPoints = finitePoints & ~outsideGrid & ~axisPoints & (psiN > 0) & (psiN < 1);

% Sample the bicubic interpolant once for all eligible points. The
% transposed array gives contourc its Z-by-R layout.
if any(interiorPoints)
    [rFine,zFine,psiNGrid] = refine_psiN(g.r,g.z,refinement,g);
    psiNGrid = psiNGrid.';
    gridSpacing = hypot(rFine(2)-rFine(1),zFine(2)-zFine(1));
end

% Locate the outboard separatrix once. The radial line also brackets the
% outboard intersection of each interior surface without using g.bdry.
if any(interiorPoints | axisPoints)
    Omp = outboard_midplane_map(g,refinement,zOMP);
end

% Each point has an independent status. A rejected point leaves other
% surfaces intact; malformed call-wide inputs are still errors above.
Geo = repmat(empty_geometry(NaN,refinement,zOMP,Reference),size(R));
for iPoint = 1:numel(rPoints)
    Geo(iPoint).psiN = psiN(iPoint);
    if ~finitePoints(iPoint)
        Geo(iPoint).errorId = 'InvalidPoint';
    elseif outsideGrid(iPoint)
        Geo(iPoint).errorId = 'OutsideGrid';
    elseif axisPoints(iPoint)
        if ~isfinite(Omp.rSep)
            Geo(iPoint).errorId = 'NoOmpSeparatrix';
        else
            Geo(iPoint).areaRZ = 0;
            Geo(iPoint).volume = 0;
            Geo(iPoint).phiTor = 0;
            Geo(iPoint).rTor = 0;
            Geo(iPoint).rSaarelma = 0;
            Geo(iPoint).dVdrTor = 0;
            Geo(iPoint).dVdrSaarelma = 0;
            Geo(iPoint).transportMetricRTor = 0;
            Geo(iPoint).transportMetricRSaarelma = 0;
            if zOMP == g.zmaxis
                Geo(iPoint).x = rPoints(iPoint)-Omp.rSep;
                Geo(iPoint).dVdx = 0;
                Geo(iPoint).dpsiNdx = 0;
                Geo(iPoint).transportMetric = 0;
            end
            Geo(iPoint).centroidR = rPoints(iPoint);
            Geo(iPoint).contourR = rPoints(iPoint);
            Geo(iPoint).contourZ = zPoints(iPoint);
            Geo(iPoint).iErr = 0;
        end
    elseif ~interiorPoints(iPoint)
        Geo(iPoint).errorId = 'OutsideCore';
    else
        Geo(iPoint) = geometry_at_point(g,rPoints(iPoint),zPoints(iPoint),psiN(iPoint),refinement,rFine,zFine,psiNGrid,gridSpacing,Omp,Reference,Fpol);
    end
end
end

function Geo = empty_geometry(psiN,refinement,zOMP,Reference)
% Use NaN for unavailable scalar geometry and empty arrays for contours.
Geo = struct('psiN',psiN,'x',NaN,'areaRZ',NaN,'volume',NaN,'dVdx',NaN,'gradXSqAvg',NaN,'dpsiNdx',NaN,'transportMetric',NaN,'centroidR',NaN,'contourR',[],'contourZ',[],'refinement',refinement,'zOMP',zOMP,'xDefinition','R_out(psiN,zOMP) - R_out(1,zOMP)','iErr',1,'errorId','');
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
Geo.b0Reference = Reference.b0;
Geo.r0Reference = Reference.r0;
Geo.toroidalFieldSign = Reference.fieldSign;
end

function Geo = geometry_at_point(g,R,Z,psiN,refinement,rFine,zFine,psiNGrid,gridSpacing,Omp,Reference,Fpol)
% Select the component nearest this point, then require a closed core loop.
% Private-flux and coil contours may share the same psiN as a core contour.
Geo = empty_geometry(psiN,refinement,Omp.z,Reference);
contours = contourc(rFine,zFine,psiNGrid,[psiN psiN]);
if isempty(contours)
    Geo.errorId = 'NoContour';
    return
end
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
    distance = point_to_contour_distance(R,Z,points);
    if distance < bestDistance
        bestDistance = distance;
        bestContour = points;
    end
end
if isempty(bestContour) || bestDistance > 0.25*gridSpacing
    Geo.errorId = 'NoMatchingContour';
    return
end
closureGap = hypot(bestContour(1,1)-bestContour(1,end),bestContour(2,1)-bestContour(2,end));
if size(bestContour,2) < 4 || closureGap > 1e-6*gridSpacing
    Geo.errorId = 'OpenContour';
    return
end
bestContour(:,end) = bestContour(:,1);
if ~inpolygon(g.rmaxis,g.zmaxis,bestContour(1,:),bestContour(2,:))
    Geo.errorId = 'WrongBranch';
    return
end

% Polygon cross products give the R-Z area and first radial moment. The
% volume of revolution is 2*pi*areaRZ*centroidR.
r0 = bestContour(1,1:end-1);
r1 = bestContour(1,2:end);
z0 = bestContour(2,1:end-1);
z1 = bestContour(2,2:end);
crossProducts = r0.*z1-r1.*z0;
signedArea = sum(crossProducts)/2;
if signedArea == 0
    Geo.errorId = 'InvalidContour';
    return
end
areaRZ = abs(signedArea);
centroidR = sum((r0+r1).*crossProducts)/(6*signedArea);
if centroidR <= 0
    Geo.errorId = 'InvalidContour';
    return
end

% On the outboard midplane x changes one-for-one with R. The derivative
% dpsiN/dx there converts the contour integrals from psiN to x.
if ~isfinite(Omp.rSep)
    Geo.errorId = 'NoOmpSeparatrix';
    return
end
[rOMP,dpsiNdx] = outboard_midplane_at_psiN(g,psiN,Omp);
if ~isfinite(rOMP)
    Geo.errorId = 'NoOmpIntersection';
    return
end
if ~isfinite(dpsiNdx) || dpsiNdx <= 0 || rOMP >= Omp.rSep
    Geo.errorId = 'InvalidOmpMap';
    return
end

% For an axisymmetric surface dS = 2*pi*R*dl. Integrate the volume
% element and metric using the same bicubic flux derivatives as calc_psiN.
segmentLength = hypot(r1-r0,z1-z0);
validSegments = segmentLength > 0;
if ~any(validSegments)
    Geo.errorId = 'InvalidContour';
    return
end
rMid = (r0(validSegments)+r1(validSegments))/2;
zMid = (z0(validSegments)+z1(validSegments))/2;
segmentLength = segmentLength(validSegments);
[dpsidr,dpsidz,ierr] = calc_psi_derivs(g,rMid,zMid,1);
gradPsiN = hypot(dpsidr,dpsidz)/abs(g.ssibry-g.ssimag);
if any(ierr | ~isfinite(gradPsiN) | gradPsiN <= 0)
    Geo.errorId = 'InvalidFluxGradient';
    return
end
integralInverse = sum(rMid.*segmentLength./gradPsiN);
integralDirect = sum(rMid.*segmentLength.*gradPsiN);
dVdx = 2*pi*dpsiNdx*integralInverse;
gradXSqAvg = integralDirect/(dpsiNdx^2*integralInverse);
transportMetric = 2*pi*integralDirect/dpsiNdx;
if ~isfinite(dVdx) || ~isfinite(gradXSqAvg) || ~isfinite(transportMetric) || dVdx <= 0 || gradXSqAvg <= 0 || transportMetric <= 0
    Geo.errorId = 'InvalidFluxGradient';
    return
end

% Integrating each enclosed cross-section includes the axis region even
% when only pedestal surfaces are requested. The local contour integral
% gives exact coordinate chain rules up to the contour quadrature accuracy.
[phiTor,inverseRArea,errorId] = enclosed_radial_integrals(g,bestContour,Fpol);
if ~isempty(errorId)
    Geo.errorId = errorId;
    return
end
rTor = sqrt(phiTor/(pi*Reference.b0));
rSaarelma = sqrt(Reference.r0*inverseRArea/pi);
integralToroidal = sum(segmentLength./(rMid.*gradPsiN));
dPhiTorDpsiN = ppval(Fpol,psiN)*integralToroidal;
drTorDx = dPhiTorDpsiN*dpsiNdx/(2*pi*Reference.b0*rTor);
drSaarelmaDx = Reference.r0*integralToroidal*dpsiNdx/(2*pi*rSaarelma);
radiusValues = [rTor,rSaarelma,dPhiTorDpsiN,drTorDx,drSaarelmaDx];
if any(~isfinite(radiusValues) | radiusValues <= 0)
    Geo.errorId = 'InvalidRadialMap';
    return
end

Geo.phiTor = phiTor;
Geo.rTor = rTor;
Geo.rSaarelma = rSaarelma;
Geo.dPhiTorDpsiN = dPhiTorDpsiN;
Geo.drTorDx = drTorDx;
Geo.drSaarelmaDx = drSaarelmaDx;
Geo.dVdrTor = dVdx/drTorDx;
Geo.dVdrSaarelma = dVdx/drSaarelmaDx;
Geo.gradRTorSqAvg = gradXSqAvg*drTorDx^2;
Geo.gradRSaarelmaSqAvg = gradXSqAvg*drSaarelmaDx^2;
Geo.transportMetricRTor = transportMetric*drTorDx;
Geo.transportMetricRSaarelma = transportMetric*drSaarelmaDx;
Geo.x = rOMP-Omp.rSep;
Geo.areaRZ = areaRZ;
Geo.volume = 2*pi*areaRZ*centroidR;
Geo.dVdx = dVdx;
Geo.gradXSqAvg = gradXSqAvg;
Geo.dpsiNdx = dpsiNdx;
Geo.transportMetric = transportMetric;
Geo.centroidR = centroidR;
Geo.contourR = bestContour(1,:);
Geo.contourZ = bestContour(2,:);
Geo.iErr = 0;
end

function [phiTor,inverseRArea,errorId] = enclosed_radial_integrals(g,contour,Fpol)
% Integrate |F(psiN)|/R and 1/R over the entire enclosed poloidal area.
% Triangulation includes concave regions without assuming radial rays from
% the magnetic axis stay inside the surface. The integrals include no 2*pi:
% phiTor is total toroidal flux, and rSaarelma^2 = R0*inverseRArea/pi.
phiTor = NaN;
inverseRArea = NaN;
errorId = 'InvalidRadialIntegral';
if any(contour(1,:) <= 0)
    return
end
Boundary = polyshape(contour(1,1:end-1),contour(2,1:end-1),'Simplify',false);
if Boundary.NumRegions ~= 1 || Boundary.NumHoles ~= 0
    return
end
Mesh = triangulation(Boundary);
vertices = Mesh.Points;
triangles = Mesh.ConnectivityList;
if isempty(triangles)
    return
end
a = vertices(triangles(:,1),:);
b = vertices(triangles(:,2),:);
c = vertices(triangles(:,3),:);
coarse = radial_triangle_integrals(g,a,b,c,Fpol);
if isempty(coarse) || any(~isfinite(coarse(:)))
    return
end

% Compare a degree-two triangle rule with the sum over four children.
% Refine only unresolved triangles; positive integrands make the summed
% local relative-error estimates a relative-error estimate for each integral.
total = [0,0];
relativeTolerance = 1e-6;
for iDepth = 1:8
    if size(a,1) > 250000
        break
    end
    ab = (a+b)/2;
    bc = (b+c)/2;
    ca = (c+a)/2;
    childA = [a;ab;ca;ab];
    childB = [ab;b;bc;bc];
    childC = [ca;bc;c;ca];
    childValues = radial_triangle_integrals(g,childA,childB,childC,Fpol);
    if any(~isfinite(childValues(:)))
        return
    end
    nTriangles = size(a,1);
    fine = childValues(1:nTriangles,:)+childValues(nTriangles+(1:nTriangles),:)+childValues(2*nTriangles+(1:nTriangles),:)+childValues(3*nTriangles+(1:nTriangles),:);
    accepted = all(abs(fine-coarse) <= relativeTolerance*max(abs(fine),realmin),2);
    total = total+sum(fine(accepted,:),1);
    if all(accepted)
        if all(isfinite(total) & total > 0)
            phiTor = total(1);
            inverseRArea = total(2);
            errorId = '';
        end
        return
    end
    unresolved = repmat(~accepted,4,1);
    a = childA(unresolved,:);
    b = childB(unresolved,:);
    c = childC(unresolved,:);
    coarse = childValues(unresolved,:);
end
errorId = 'RadialIntegralNotConverged';
end

function values = radial_triangle_integrals(g,a,b,c,Fpol)
% Three equal-weight barycentric points integrate quadratic polynomials.
% Each row contains integrals of |F|/R and 1/R over one physical triangle.
area = abs((b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(c(:,1)-a(:,1)).*(b(:,2)-a(:,2)))/2;
values = zeros(size(a,1),2);
barycentric = [2/3,1/6,1/6;1/6,2/3,1/6;1/6,1/6,2/3];
for iPoint = 1:3
    points = barycentric(iPoint,1)*a+barycentric(iPoint,2)*b+barycentric(iPoint,3)*c;
    [psiN,~,ierr] = calc_psiN(g,points(:,1).',points(:,2).',1);
    if any(ierr | ~isfinite(psiN)) || any(points(:,1) <= 0)
        values(:) = NaN;
        return
    end
    % Endpoint values of F handle flux overshoots from representing a
    % smooth contour by polygon edges and interpolating the gfile near its axis.
    fMagnitude = ppval(Fpol,max(0,min(1,psiN))).';
    values = values+[fMagnitude./points(:,1),1./points(:,1)]/3;
end
values = values.*area;
end

function Omp = outboard_midplane_map(g,refinement,zOMP)
% Bracket psiN = 1 on the line from the axis radius toward the outboard edge.
% The first crossing excludes unrelated flux structure farther outside.
rLine = linspace(g.rmaxis,g.r(end),max(256,refinement*numel(g.r)));
psiNLine = calc_psiN(g,rLine,zOMP*ones(size(rLine)),1);
iSep = find(psiNLine(1:end-1) < 1 & psiNLine(2:end) >= 1,1,'first');
Omp = struct('z',zOMP,'rLine',rLine,'psiNLine',psiNLine,'iSep',0,'rSep',NaN);
if ~isempty(iSep)
    Omp.iSep = iSep;
    Omp.rSep = fzero(@(rValue) calc_psiN(g,rValue,zOMP,1)-1,rLine(iSep:iSep+1));
end
end

function [rOMP,dpsiNdx] = outboard_midplane_at_psiN(g,psiN,Omp)
% Select the outermost outward crossing inside the OMP separatrix, then
% evaluate dpsiN/dR from the bicubic interpolant at that intersection.
rOMP = NaN;
dpsiNdx = NaN;
if Omp.iSep == 0
    return
end
iCross = find(Omp.psiNLine(1:Omp.iSep) <= psiN & Omp.psiNLine(2:Omp.iSep+1) >= psiN,1,'last');
if isempty(iCross)
    return
end
rOMP = fzero(@(rValue) calc_psiN(g,rValue,Omp.z,1)-psiN,Omp.rLine(iCross:iCross+1));
[dpsidr,~,ierr] = calc_psi_derivs(g,rOMP,Omp.z,1);
if ~ierr
    % calc_psi_derivs includes ip_sign; calc_psiN applies it a second time.
    dpsiNdx = g.ip_sign*dpsidr/(g.ssibry-g.ssimag);
end
end

function distance = point_to_contour_distance(R,Z,points)
% Find the shortest distance to any line segment, including its endpoints.
startPoints = points(:,1:end-1);
segments = points(:,2:end)-startPoints;
lengthSquared = sum(segments.^2,1);
lengthSquared(lengthSquared == 0) = 1;
fraction = ((R-startPoints(1,:)).*segments(1,:)+(Z-startPoints(2,:)).*segments(2,:))./lengthSquared;
fraction = max(0,min(1,fraction));
closest = startPoints+segments.*fraction;
distance = min(hypot(R-closest(1,:),Z-closest(2,:)));
end
