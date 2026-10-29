function test_compute_flux_surface_geometry
% Check geometry and contour selection on analytic, disconnected flux grids.

testDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(testDirectory,'..','..','bfield_library_jdl','EFIT_routines'));

% A quadratic flux field has elliptical level sets with known area and
% volume. Leave bdry empty to ensure the routine does not require it.
rAxis = 2;
a = 0.5;
b = 0.7;
g = struct('mw',65,'mh',65,'xdim',2,'zdim',2,'rgrid1',1,'zmid',0,'rmaxis',rAxis,'zmaxis',0,'ssimag',0,'ssibry',1,'cpasma',[],'bdry',[],'fpol',ones(1,65));
[rGrid,zGrid] = ndgrid(linspace(1,3,g.mw),linspace(-1,1,g.mh));
corePsi = ((rGrid-rAxis)/a).^2+(zGrid/b).^2;
g.psirz = corePsi;
g = postprocess_gfile(g);
Axis = compute_flux_surface_geometry(g,rAxis,0);
assert(Axis.areaRZ == 0 && Axis.volume == 0 && Axis.refinement == 4 && Axis.iErr == 0);
assert(abs(Axis.x+a) < 1e-8 && Axis.dVdx == 0 && Axis.dpsiNdx == 0);
assert(isnan(Axis.gradXSqAvg));

for psiN = [0.25,0.8,0.98]
    Geometry = compute_flux_surface_geometry(g,rAxis+a*sqrt(psiN),0);
    GeometryFine = compute_flux_surface_geometry(g,rAxis+a*sqrt(psiN),0,6);
    expectedArea = pi*a*b*psiN;
    expectedVolume = 2*pi*rAxis*expectedArea;
    expectedX = a*(sqrt(psiN)-1);
    expectedDVdx = 4*pi^2*rAxis*b*sqrt(psiN);
    expectedGradXSqAvg = (1+(a/b)^2)/2;
    assert(abs(Geometry.areaRZ/expectedArea-1) < 0.01);
    assert(abs(Geometry.volume/expectedVolume-1) < 0.01);
    assert(abs(Geometry.x-expectedX) < 1e-6);
    assert(abs(Geometry.dVdx/expectedDVdx-1) < 0.02);
    assert(abs(Geometry.gradXSqAvg/expectedGradXSqAvg-1) < 0.02);
    assert(abs(Geometry.dpsiNdx-2*sqrt(psiN)/a) < 1e-5);
    assert(abs(Geometry.centroidR-rAxis) < 0.01);
    assert(abs(GeometryFine.areaRZ/Geometry.areaRZ-1) < 0.01);
    assert(GeometryFine.refinement == 6 && Geometry.iErr == 0 && GeometryFine.iErr == 0);
end

% The precomputed table should reproduce smooth analytic geometry while
% resolving x at Z = 0 independently of the line used to select levels.
gTable = postprocess_gfile_geometry(g,false,4,[0.1,0.25,0.5,0.75,0.95]);
GeoLookup = lookup_flux_surface_geometry(gTable,rAxis+a*sqrt(0.6),0);
assert(GeoLookup.iErr == 0 && GeoLookup.zOMP == 0 && isempty(GeoLookup.contourR));
assert(abs(GeoLookup.x-a*(sqrt(0.6)-1)) < 0.01);
assert(abs(GeoLookup.volume-2*pi^2*rAxis*a*b*0.6) < 0.02);
assert(abs(GeoLookup.gradXSqAvg-(1+(a/b)^2)/2) < 0.02);
assert(strcmp(gTable.fluxSurfaceGeometry.xDefinition,'R_out(psiN,Z=0) - R_out(1,Z=0)'));
assert(gTable.fluxSurfaceGeometry.zSampling == g.zmaxis);
% A sparse table should retain the curved x(psiN) relation between nodes.
gSparse = postprocess_gfile_geometry(g,false,4,[0.1,0.3,0.6,0.9]);
psiNSparse = 0.45;
GeoSparse = lookup_flux_surface_geometry(gSparse,rAxis+a*sqrt(psiNSparse),0);
assert(GeoSparse.iErr == 0 && abs(GeoSparse.x-a*(sqrt(psiNSparse)-1)) < 0.001);
GeoLookupMixed = lookup_flux_surface_geometry(gTable,[rAxis+a*sqrt(0.6),rAxis+a*sqrt(0.99)],0);
assert(GeoLookupMixed(1).iErr == 0 && strcmp(GeoLookupMixed(2).errorId,'OutsideTable'));
assert_status(lookup_flux_surface_geometry(gTable,rAxis+a*sqrt(0.99),0),'OutsideTable');

% The optional geometry field should be written into, then reused from,
% the existing version-10 gfile cache without changing its base version.
cacheBase = tempname;
cacheFile = [cacheBase,'.mat'];
gCache = g;
gCache.filename = cacheBase;
gCache.version = 10;
Cache = struct('g',gCache);
save(cacheFile,'-struct','Cache');
gCache = postprocess_gfile_geometry(gCache,false,4,[0.1,0.25,0.5,0.75,0.95]);
Cache = load(cacheFile,'g');
assert(Cache.g.version == 10 && isfield(Cache.g,'fluxSurfaceGeometry'));
Cache.g.fluxSurfaceGeometry.reuseMarker = true;
save(cacheFile,'-struct','Cache');
gCache = postprocess_gfile_geometry(Cache.g,false,4,[0.1,0.25,0.5,0.75,0.95]);
assert(gCache.fluxSurfaceGeometry.reuseMarker);
delete(cacheFile);

gShift = g;
gShift.zmaxis = 0.2;
gShift.psirz = ((rGrid-rAxis)/a).^2+((zGrid-gShift.zmaxis)/b).^2;
gShift = postprocess_gfile(gShift);
gShift = postprocess_gfile_geometry(gShift,false,4,[0.1,0.25,0.5,0.75,0.9]);
GeoShift = lookup_flux_surface_geometry(gShift,rAxis+a*sqrt(0.8),gShift.zmaxis);
expectedShiftX = a*(sqrt(0.8-(gShift.zmaxis/b)^2)-sqrt(1-(gShift.zmaxis/b)^2));
assert(GeoShift.iErr == 0 && abs(GeoShift.x-expectedShiftX) < 0.01);

% Matching R and Z arrays return a Geo struct array of the same shape.
% A scalar Z expands across R, and each result matches a scalar call.
rPoints = [rAxis,rAxis+a*sqrt(0.25);rAxis+a*sqrt(0.8),rAxis+a*sqrt(0.98)];
GeoArray = compute_flux_surface_geometry(g,rPoints,0);
assert(isequal(size(GeoArray),size(rPoints)));
for iPoint = 1:numel(rPoints)
    GeoScalar = compute_flux_surface_geometry(g,rPoints(iPoint),0,4);
    assert(GeoArray(iPoint).areaRZ == GeoScalar.areaRZ);
    assert(isequal(GeoArray(iPoint).contourR,GeoScalar.contourR));
    assert(GeoArray(iPoint).refinement == 4 && GeoArray(iPoint).iErr == 0);
end
GeoArrayFine = compute_flux_surface_geometry(g,rPoints,zeros(size(rPoints)),6);
assert(isequal(size(GeoArrayFine),size(rPoints)));
assert(all([GeoArrayFine.refinement] == 6) && all([GeoArrayFine.iErr] == 0));
assert_error(@() compute_flux_surface_geometry(g,[rAxis,rAxis],zeros(2,1)),'compute_flux_surface_geometry:InvalidPoint');
GeoMixed = compute_flux_surface_geometry(g,[rAxis,rAxis+0.55,4,NaN,rAxis+a*sqrt(0.8)],0);
assert(GeoMixed(1).iErr == 0 && GeoMixed(1).areaRZ == 0);
assert_status(GeoMixed(2),'OutsideCore');
assert(GeoMixed(2).psiN > 1);
assert_status(GeoMixed(3),'OutsideGrid');
assert_status(GeoMixed(4),'InvalidPoint');
assert(isnan(GeoMixed(3).psiN) && isnan(GeoMixed(4).psiN));
assert(GeoMixed(5).iErr == 0 && GeoMixed(5).areaRZ > 0);

% The 65-by-65 grid needs factor 4 to reach 256 points. A taller grid
% needs factor 5 to keep the refined Z spacing below 5 cm.
GeometryAuto = compute_flux_surface_geometry(g,rAxis+a*sqrt(0.8),0);
GeometryFour = compute_flux_surface_geometry(g,rAxis+a*sqrt(0.8),0,4);
assert(GeometryAuto.refinement == 4 && GeometryFour.refinement == 4);
assert(isequal(GeometryAuto.contourR,GeometryFour.contourR));
gTall = g;
gTall.zdim = 14;
[~,zTall] = ndgrid(g.r,linspace(-7,7,g.mh));
gTall.psirz = ((rGrid-rAxis)/a).^2+(zTall/b).^2;
gTall = postprocess_gfile(gTall);
GeometryAuto = compute_flux_surface_geometry(gTall,rAxis+a*sqrt(0.8),0);
GeometryFive = compute_flux_surface_geometry(gTall,rAxis+a*sqrt(0.8),0,5);
assert(GeometryAuto.refinement == 5 && GeometryFive.refinement == 5);
assert(isequal(GeometryAuto.contourR,GeometryFive.contourR));

% A flux contour clipped by the computational box must not be treated as
% closed, even when its normalized flux value is below one.
gOpen = g;
gOpen.ssibry = 10;
assert_status(compute_flux_surface_geometry(gOpen,2.8,0),'OpenContour');

% A second low-flux island has the same psiN as the core but does not
% surround the magnetic axis. The query must not inherit the core geometry.
g.psirz = min(corePsi,((rGrid-2.8)/0.15).^2+(zGrid/0.15).^2);
g = postprocess_gfile(g);
assert_status(compute_flux_surface_geometry(g,2.8+0.15*sqrt(0.25),0),'WrongBranch');
assert_status(compute_flux_surface_geometry(g,rAxis+0.55,0),'OutsideCore');
gIslandTable = postprocess_gfile_geometry(g,false,4,[0.1,0.25,0.5,0.75,0.9]);
assert_status(lookup_flux_surface_geometry(gIslandTable,2.8+0.15*sqrt(0.25),0),'WrongBranch');

fprintf('test_compute_flux_surface_geometry passed\n');
end

function assert_status(Geo,expectedId)
% Check that a rejected point has a useful status and no geometry.
assert(Geo.iErr == 1 && strcmp(Geo.errorId,expectedId));
assert(isnan(Geo.areaRZ) && isnan(Geo.volume) && isnan(Geo.centroidR));
assert(isnan(Geo.x) && isnan(Geo.dVdx) && isnan(Geo.gradXSqAvg) && isnan(Geo.dpsiNdx));
assert(isempty(Geo.contourR) && isempty(Geo.contourZ));
end

function assert_error(callFunction,expectedIdentifier)
% Check the reason for rejection without depending on warning text.
try
    callFunction();
catch ME
    assert(strcmp(ME.identifier,expectedIdentifier),ME.message);
    return
end
error('Expected error %s was not raised',expectedIdentifier)
end
