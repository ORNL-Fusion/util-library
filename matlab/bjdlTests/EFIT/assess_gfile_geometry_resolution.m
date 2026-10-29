function Results = assess_gfile_geometry_resolution(gfileName,baseRefinement,fineRefinement)
% Compare geometry-table interpolation and contour-grid resolution.
% Results = assess_gfile_geometry_resolution(gfileName,baseRefinement,fineRefinement)
% The default gfileName is the MAST-U case under C:\Work. Empty refinement
% inputs select the production default and twice its contour factor.
% This diagnostic plots errors at flux levels halfway between table nodes.

if nargin < 1 || isempty(gfileName)
    gfileName = 'C:\Work\SOLPS\runs\mastu\cd-51685-400ms\baserun\gfile';
end
if nargin < 2
    baseRefinement = [];
end
if nargin < 3
    fineRefinement = [];
end

% Load the equilibrium and build a fresh table without saving this
% diagnostic's optional geometry into the normal gfile cache.
testDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(testDirectory,'..','..','bfield_library_jdl','EFIT_routines'));
addpath(fullfile(testDirectory,'..','..','bfield_library_jdl','util'));
g = readg_g3d(gfileName);
if isempty(g)
    error('assess_gfile_geometry_resolution:NoGfile','Could not read %s',gfileName)
end
if isfield(g,'fluxSurfaceGeometry')
    g = rmfield(g,'fluxSurfaceGeometry');
end
g.filename = '';
g = postprocess_gfile_geometry(g,false,baseRefinement);
Table = g.fluxSurfaceGeometry;
baseRefinement = Table.refinement;
if isempty(fineRefinement)
    fineRefinement = 2*baseRefinement;
end
if ~isscalar(fineRefinement) || ~isnumeric(fineRefinement) || ~isfinite(fineRefinement) || fineRefinement <= baseRefinement || fineRefinement ~= floor(fineRefinement)
    error('assess_gfile_geometry_resolution:InvalidRefinement','fineRefinement must be an integer greater than %d',baseRefinement)
end

% Locate independent test points halfway between retained table levels.
% The axis-height line chooses a point; each direct call uses its full
% closed contour and defines x on Z = 0, just like the table builder.
psiNMid = (Table.psiN(1:end-1)+Table.psiN(2:end))/2;
rMid = radii_at_flux_levels(g,psiNMid);
sampled = isfinite(rMid);
if ~all(sampled)
    warning('assess_gfile_geometry_resolution:NoMidpoint','Could not locate %d midpoint levels on the axis-height line',nnz(~sampled))
end
psiNMid = psiNMid(sampled);
rMid = rMid(sampled);
if isempty(rMid)
    error('assess_gfile_geometry_resolution:NoMidpoint','No midpoint levels could be located')
end

% Compare the PCHIP lookup with direct contours on the production grid,
% then repeat the direct calculation on a finer contour grid.
DirectBase = compute_flux_surface_geometry(g,rMid,g.zmaxis,baseRefinement,0);
DirectFine = compute_flux_surface_geometry(g,rMid,g.zmaxis,fineRefinement,0);
Lookup = lookup_flux_surface_geometry(g,rMid,g.zmaxis);
valid = [DirectBase.iErr] == 0 & [DirectFine.iErr] == 0 & [Lookup.iErr] == 0;
fprintf('\nGeometry resolution for %s\n',gfileName)
fprintf('Requested/retained levels: %d/%d\n',Table.nRequested,Table.nComputed)
skippedLevels = Table.requestedPsiNLevels(arrayfun(@(level) all(abs(Table.psiN-level) > 1e-8),Table.requestedPsiNLevels));
if ~isempty(skippedLevels)
    fprintf('Skipped psiN levels: %s\n',mat2str(skippedLevels,7))
end
fprintf('Base/fine contour refinement: %d/%d\n',baseRefinement,fineRefinement)
fprintf('Valid midpoint comparisons: %d/%d\n',nnz(valid),numel(psiNMid))
for iPoint = find(~valid)
    fprintf('  psiN = %.8f: direct %s, fine %s, lookup %s\n',psiNMid(iPoint),DirectBase(iPoint).errorId,DirectFine(iPoint).errorId,Lookup(iPoint).errorId)
end
if ~any(valid)
    error('assess_gfile_geometry_resolution:NoValidMidpoints','No midpoint has valid geometry at both contour refinements')
end

% Scale discrepancies by the direct value with a 1%% profile-scale floor.
% This keeps errors meaningful for x and volume close to zero.
fieldNames = {'x','areaRZ','volume','dVdx','gradXSqAvg','transportMetric'};
Results = struct();
Results.gfileName = gfileName;
Results.psiN = psiNMid;
Results.r = rMid;
Results.tablePsiN = Table.psiN;
Results.skippedLevels = skippedLevels;
Results.baseRefinement = baseRefinement;
Results.fineRefinement = fineRefinement;
Results.valid = valid;
Results.lookupError = struct();
Results.refinementError = struct();
fprintf('Scaled maximum discrepancies (denominator: max(|direct|, 1%% of profile scale))\n')
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    baseValues = [DirectBase.(fieldName)];
    fineValues = [DirectFine.(fieldName)];
    lookupValues = [Lookup.(fieldName)];
    fieldValid = valid & isfinite(baseValues) & isfinite(fineValues) & isfinite(lookupValues);
    profileScale = max(abs(Table.(fieldName)));
    denominator = max(abs(baseValues),0.01*max(profileScale,eps));
    lookupError = nan(size(psiNMid));
    refinementError = nan(size(psiNMid));
    lookupError(fieldValid) = abs(lookupValues(fieldValid)-baseValues(fieldValid))./denominator(fieldValid);
    refinementError(fieldValid) = abs(fineValues(fieldValid)-baseValues(fieldValid))./denominator(fieldValid);
    Results.lookupError.(fieldName) = lookupError;
    Results.refinementError.(fieldName) = refinementError;
    if any(fieldValid)
        fprintf('  %-16s lookup %.3g, contour refinement %.3g\n',fieldName,max(lookupError(fieldValid)),max(refinementError(fieldValid)))
    else
        fprintf('  %-16s no valid comparisons\n',fieldName)
    end
end

% Plot both error sources against flux to expose edge or axis trends.
figure('Color','w','Name','Gfile geometry resolution errors');
tiledlayout(2,3);
for iField = 1:numel(fieldNames)
    fieldName = fieldNames{iField};
    nexttile;
    semilogy(psiNMid,max(Results.lookupError.(fieldName),eps),'b.-',psiNMid,max(Results.refinementError.(fieldName),eps),'r.-');
    grid on; xlim([0,1]); xlabel('\psi_N'); ylabel('Scaled discrepancy'); title(fieldName,'Interpreter','none');
    if iField == 1
        legend('Lookup versus direct','Fine versus base','Location','best');
    end
end
sgtitle(sprintf('Geometry resolution: contour factors %d and %d',baseRefinement,fineRefinement));

% Show the actual sampling density and the edge behavior of dV/dx.
figure('Color','w','Name','Gfile geometry sampling');
subplot(2,1,1);
plot(Table.psiN,Table.dVdx,'k.-',psiNMid(valid),[DirectFine(valid).dVdx],'ro');
grid on; xlim([0,1]); xlabel('\psi_N'); ylabel('dV/dx (m^2)');
legend('Stored table','Fine-grid midpoints','Location','best');
subplot(2,1,2);
semilogy(Table.psiN(1:end-1),diff(Table.psiN),'k.-');
grid on; xlim([0,1]); xlabel('\psi_N'); ylabel('Table spacing in \psi_N');
end

function radii = radii_at_flux_levels(g,levels)
% Bracket the first outward crossing of each level before the separatrix.
rLine = linspace(g.rmaxis,g.r(end),max(512,2*numel(g.r)));
psiNLine = calc_psiN(g,rLine,g.zmaxis*ones(size(rLine)),1);
iSep = find(psiNLine(1:end-1) < 1 & psiNLine(2:end) >= 1,1,'first');
if isempty(iSep)
    error('assess_gfile_geometry_resolution:NoSeparatrix','No outboard separatrix crossing was found at the axis height')
end
radii = nan(size(levels));
for iLevel = 1:numel(levels)
    iCross = find(psiNLine(1:iSep) <= levels(iLevel) & psiNLine(2:iSep+1) >= levels(iLevel),1,'first');
    if ~isempty(iCross)
        radii(iLevel) = fzero(@(rValue) calc_psiN(g,rValue,g.zmaxis,1)-levels(iLevel),rLine(iCross:iCross+1));
    end
end
end
