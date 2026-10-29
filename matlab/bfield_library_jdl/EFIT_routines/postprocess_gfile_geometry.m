function g = postprocess_gfile_geometry(g,makePlots,refinement,psiNLevels)
% Add a reusable table of closed-core flux-surface geometry to g.
% Flux levels are selected on the line Z = g.zmaxis, but x is defined by
% the outboard intersections at Z = 0: R_out(psiN,0)-R_out(1,0).
% Usage: g = postprocess_gfile_geometry(g,makePlots,refinement,psiNLevels)
% makePlots is optional (default false); refinement and psiNLevels are
% optional overrides. A valid table is reused, including when readg_g3d
% loaded it from the gfile .mat cache. Newly computed tables update that cache.

if nargin < 2 || isempty(makePlots)
    makePlots = false;
end
if nargin < 3
    refinement = [];
end
if nargin < 4 || isempty(psiNLevels)
    % Resolve the smooth interior sparsely and approach the separatrix
    % geometrically. Unresolved levels are discarded after contour checks.
    psiNLevels = unique([linspace(0.02,0.9,24),1-logspace(-1,-4,20)]);
end
if ~isscalar(makePlots) || (~islogical(makePlots) && ~isnumeric(makePlots)) || ~isreal(makePlots) || ~isfinite(makePlots) || ~ismember(makePlots,[0,1])
    error('postprocess_gfile_geometry:InvalidPlotFlag','makePlots must be true or false')
end
if ~isnumeric(psiNLevels) || ~isvector(psiNLevels) || numel(psiNLevels) < 2 || ~isreal(psiNLevels) || any(~isfinite(psiNLevels)) || any(psiNLevels <= 0 | psiNLevels >= 1)
    error('postprocess_gfile_geometry:InvalidLevels','psiNLevels must contain at least two finite interior levels')
end
psiNLevels = unique(psiNLevels(:).');
if ~isstruct(g) || ~all(isfield(g,{'r','z','rmaxis','zmaxis','bicub_coeffs_inv'})) || g.z(1) > 0 || g.z(end) < 0
    error('postprocess_gfile_geometry:InvalidGfile','g must be postprocessed and its grid must include Z = 0')
end
if isfield(g,'fluxSurfaceGeometry') && valid_geometry_table(g.fluxSurfaceGeometry,g,refinement,psiNLevels)
    if makePlots
        plot_gfile_geometry_table(g);
    end
    return
end

% Find points at the requested levels on the outboard line through the
% magnetic axis. That line chooses levels only; the geometry routine uses
% whole contours and the separate Z = 0 line to define x.
rLine = linspace(g.rmaxis,g.r(end),max(512,2*numel(g.r)));
psiNLine = calc_psiN(g,rLine,g.zmaxis*ones(size(rLine)),1);
iSep = find(psiNLine(1:end-1) < 1 & psiNLine(2:end) >= 1,1,'first');
if isempty(iSep)
    error('postprocess_gfile_geometry:NoSeparatrix','No outboard separatrix crossing was found at the axis height')
end
rSamples = nan(size(psiNLevels));
for iLevel = 1:numel(psiNLevels)
    iCross = find(psiNLine(1:iSep) <= psiNLevels(iLevel) & psiNLine(2:iSep+1) >= psiNLevels(iLevel),1,'first');
    if ~isempty(iCross)
        rSamples(iLevel) = fzero(@(rValue) calc_psiN(g,rValue,g.zmaxis,1)-psiNLevels(iLevel),rLine(iCross:iCross+1));
    end
end
sampled = isfinite(rSamples);
if nnz(sampled) < 2
    error('postprocess_gfile_geometry:NoSamples','Too few requested core flux levels intersect the axis-height line')
end

% A vector call shares the refined flux grid across all levels. Retain
% only closed core contours that also have a valid Z = 0 OMP intersection.
GeoSamples = compute_flux_surface_geometry(g,rSamples(sampled),g.zmaxis,refinement,0);
good = [GeoSamples.iErr] == 0;
if nnz(good) < 2
    error('postprocess_gfile_geometry:NoSamples','Too few levels have valid closed-core geometry at Z = 0')
end
GeoSamples = GeoSamples(good);
if any(diff([GeoSamples.x]) <= 0) || any(diff([GeoSamples.volume]) <= 0)
    error('postprocess_gfile_geometry:NonmonotonicTable','x or volume is not monotonic; adjust refinement or the sampled levels')
end
if numel(GeoSamples) < numel(psiNLevels)
    warning('postprocess_gfile_geometry:SkippedLevels','Retained %d of %d requested levels with closed contours and a valid Z = 0 OMP intersection',numel(GeoSamples),numel(psiNLevels))
end

% Keep contour polygons for the core-branch check during fast lookup.
% The product V'*<|grad x|^2> is useful directly in transport equations.
Table = struct();
Table.version = 1;
Table.requestedRefinement = refinement;
Table.requestedPsiNLevels = psiNLevels;
Table.psiN = [GeoSamples.psiN];
Table.x = [GeoSamples.x];
Table.areaRZ = [GeoSamples.areaRZ];
Table.volume = [GeoSamples.volume];
Table.dVdx = [GeoSamples.dVdx];
Table.dVdxFromVolume = gradient(Table.volume,Table.x);
Table.gradXSqAvg = [GeoSamples.gradXSqAvg];
Table.dpsiNdx = [GeoSamples.dpsiNdx];
Table.transportMetric = [GeoSamples.transportMetric];
Table.centroidR = [GeoSamples.centroidR];
Table.contourR = {GeoSamples.contourR};
Table.contourZ = {GeoSamples.contourZ};
Table.refinement = GeoSamples(1).refinement;
Table.zSampling = g.zmaxis;
Table.zOMP = 0;
Table.xDefinition = 'R_out(psiN,Z=0) - R_out(1,Z=0)';
Table.nRequested = numel(psiNLevels);
Table.nComputed = numel(Table.psiN);
g.fluxSurfaceGeometry = Table;
save_geometry_to_gfile_cache(g);

if makePlots
    plot_gfile_geometry_table(g);
end
end

function valid = valid_geometry_table(Table,g,refinement,psiNLevels)
% Keep the base gfile cache version independent of optional geometry data.
fields = {'version','requestedRefinement','requestedPsiNLevels','psiN','x','areaRZ','volume','dVdx','gradXSqAvg','dpsiNdx','transportMetric','centroidR','contourR','contourZ','zSampling','zOMP'};
valid = isstruct(Table) && all(isfield(Table,fields));
if valid
    valid = Table.version == 1 && isequaln(Table.requestedRefinement,refinement) && isequaln(Table.requestedPsiNLevels,psiNLevels) && Table.zSampling == g.zmaxis && Table.zOMP == 0 && numel(Table.psiN) >= 2 && all(diff(Table.psiN) > 0);
end
end

function save_geometry_to_gfile_cache(g)
% Add the table only when this g still represents the cached raw equilibrium.
if ~isfield(g,'filename') || isempty(g.filename)
    return
end
cacheFile = [g.filename,'.mat'];
if exist(cacheFile,'file') ~= 2
    return
end
if exist(g.filename,'file') == 2
    rawInfo = dir(g.filename);
    cacheInfo = dir(cacheFile);
    if rawInfo.datenum > cacheInfo.datenum
        return
    end
end
try
    Cache = load(cacheFile,'g');
    if ~isfield(Cache,'g') || ~isfield(Cache.g,'version') || Cache.g.version ~= 10
        return
    end
    sourceFields = {'psirz','ssimag','ssibry','rmaxis','zmaxis','r','z','ip_sign'};
    if ~all(isfield(Cache.g,sourceFields)) || ~all(isfield(g,sourceFields))
        return
    end
    for iField = 1:numel(sourceFields)
        fieldName = sourceFields{iField};
        if ~isequaln(Cache.g.(fieldName),g.(fieldName))
            warning('postprocess_gfile_geometry:CacheNotUpdated','Geometry was computed on a modified gfile; the original cache was left unchanged')
            return
        end
    end
    Cache.g.fluxSurfaceGeometry = g.fluxSurfaceGeometry;
    save(cacheFile,'-struct','Cache');
catch ME
    warning('postprocess_gfile_geometry:CacheNotUpdated','Could not update geometry in %s: %s',cacheFile,ME.message)
end
end

function plot_gfile_geometry_table(g)
% Show the radial profiles and sampled contours with the limiter outline.
Table = g.fluxSurfaceGeometry;
figure('Color','w');
subplot(2,4,1); plot(Table.psiN,Table.x,'.-'); grid on; xlabel('\psi_N'); ylabel('x (m)');
subplot(2,4,2); plot(Table.x,Table.areaRZ,'.-'); grid on; xlabel('x (m)'); ylabel('R-Z area (m^2)');
subplot(2,4,3); plot(Table.x,Table.volume,'.-'); grid on; xlabel('x (m)'); ylabel('V (m^3)');
subplot(2,4,4); plot(Table.x,Table.dVdx,'.-',Table.x,Table.dVdxFromVolume,'--'); grid on; xlabel('x (m)'); ylabel('dV/dx (m^2)'); legend('Contour integral','Slope of V(x)');
subplot(2,4,5); plot(Table.x,Table.gradXSqAvg,'.-'); grid on; xlabel('x (m)'); ylabel('<|grad x|^2>');
subplot(2,4,6); plot(Table.x,Table.transportMetric,'.-'); grid on; xlabel('x (m)'); ylabel('V'' <|grad x|^2> (m^2)');
subplot(2,4,7); plot(Table.x,Table.dpsiNdx,'.-'); grid on; xlabel('x (m)'); ylabel('d\psi_N/dx (m^{-1})');
subplot(2,4,8); plot(Table.x,Table.centroidR,'.-'); grid on; xlabel('x (m)'); ylabel('R centroid (m)');

figure('Color','w'); hold on;
colors = parula(Table.nComputed);
for iLevel = 1:Table.nComputed
    plot(Table.contourR{iLevel},Table.contourZ{iLevel},'Color',colors(iLevel,:));
end
if isfield(g,'lim') && ~isempty(g.lim) && size(g.lim,1) == 2
    valid = isfinite(g.lim(1,:)) & isfinite(g.lim(2,:)) & g.lim(1,:) > 0;
    plot(g.lim(1,valid),g.lim(2,valid),'k-','LineWidth',1.5);
end
plot(g.rmaxis,g.zmaxis,'ko','MarkerFaceColor','k');
plot([g.r(1),g.r(end)],[0,0],'k:');
axis equal; grid on; xlabel('R (m)'); ylabel('Z (m)');
title('Sampled closed-core contours and limiter');
end
