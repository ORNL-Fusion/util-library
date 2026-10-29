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

requiredFields = {'r','z','psirz','rmaxis','zmaxis','ssimag','ssibry','ip_sign','bicub_coeffs_inv'};
if ~isstruct(g) || ~all(isfield(g,requiredFields))
    error('compute_flux_surface_geometry:InvalidGfile','g must be a postprocessed gfile')
end
if g.ssimag == g.ssibry
    error('compute_flux_surface_geometry:InvalidGfile','Axis and boundary flux must differ')
end
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
Geo = repmat(empty_geometry(NaN,refinement,zOMP),size(R));
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
        Geo(iPoint) = geometry_at_point(g,rPoints(iPoint),zPoints(iPoint),psiN(iPoint),refinement,rFine,zFine,psiNGrid,gridSpacing,Omp);
    end
end
end

function Geo = empty_geometry(psiN,refinement,zOMP)
% Use NaN for unavailable scalar geometry and empty arrays for contours.
Geo = struct('psiN',psiN,'x',NaN,'areaRZ',NaN,'volume',NaN,'dVdx',NaN,'gradXSqAvg',NaN,'dpsiNdx',NaN,'transportMetric',NaN,'centroidR',NaN,'contourR',[],'contourZ',[],'refinement',refinement,'zOMP',zOMP,'xDefinition','R_out(psiN,zOMP) - R_out(1,zOMP)','iErr',1,'errorId','');
end

function Geo = geometry_at_point(g,R,Z,psiN,refinement,rFine,zFine,psiNGrid,gridSpacing,Omp)
% Select the component nearest this point, then require a closed core loop.
% Private-flux and coil contours may share the same psiN as a core contour.
Geo = empty_geometry(psiN,refinement,Omp.z);
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
