function Results = assess_flux_surface_coordinates(gfileName,makePlots)
% Compare radial labels, metrics and derivatives on pedestal flux surfaces.
% Results = assess_flux_surface_coordinates(gfileName,makePlots)
% makePlots defaults to true. Read the raw equilibrium without changing its
% .mat cache. Results includes direct, cached and finer-contour calculations,
% finite-difference discrepancies, transformation residuals and plot handles.
% Relative errors use the direct calculation as the reference. These are
% numerical consistency checks, not a validation of a transport model.
if nargin < 2
    makePlots = true;
end
testDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(testDirectory,'..','..','bfield_library_jdl','EFIT_routines'));
g = postprocess_gfile(readg_g3d_simple(gfileName));
g.filename = '';
fprintf('Building geometry table for %s\n',gfileName);
g = postprocess_gfile_geometry(g,false);
Table = g.fluxSurfaceGeometry;

% Midpoints test interpolation independently of stored nodes; fixed levels
% provide reproducible comparisons across the pedestal and towards the X-point.
psiN = [0.85,0.90,0.95,0.98,0.99,0.995,0.999];
midpoints = (Table.psiN(1:end-1)+Table.psiN(2:end))/2;
queryPsi = unique([psiN,midpoints(midpoints >= 0.85 & midpoints <= 0.999)]);
assert(all(queryPsi > Table.psiN(1) & queryPsi < Table.psiN(end)),'The retained table does not cover the diagnostic interval');
fprintf('Computing %d independent direct surfaces\n',numel(queryPsi));
rQuery = radii_at_flux_levels(g,queryPsi);
Direct = compute_flux_surface_geometry(g,rQuery,g.zmaxis,Table.refinement,0);
Cached = lookup_flux_surface_geometry(g,rQuery,g.zmaxis);
assert_valid(Direct,'direct');
assert_valid(Cached,'cached');
Labels = lookup_flux_surface_labels(g,queryPsi);
Spatial = calc_flux_surface_labels(g,rQuery,g.zmaxis);
assert_valid(Labels,'flux labels');
assert_valid(Spatial,'spatial labels');
assert(max(abs([Labels.rTor]-[Spatial.rTor])) < 1e-8,'Spatial and flux label queries disagree');

% Refining the contour grid separates contour discretization from table
% interpolation. Finite differences of direct enclosed quantities provide
% a check independent of the formulas for boundary-integral derivatives.
rFixed = radii_at_flux_levels(g,psiN);
Base = compute_flux_surface_geometry(g,rFixed,g.zmaxis,Table.refinement,0);
Fine = compute_flux_surface_geometry(g,rFixed,g.zmaxis,2*Table.refinement,0);
assert_valid(Base,'base');
assert_valid(Fine,'fine');
fprintf('Computing direct finite differences at two step sizes\n');
step = min(2e-4,0.05*(1-psiN));
DirectFD = struct();
for iStep = 1:2
    h = step/2^(iStep-1);
    rPlus = radii_at_flux_levels(g,psiN+h);
    rMinus = radii_at_flux_levels(g,psiN-h);
    Plus = compute_flux_surface_geometry(g,rPlus,g.zmaxis,2*Table.refinement,0);
    Minus = compute_flux_surface_geometry(g,rMinus,g.zmaxis,2*Table.refinement,0);
    assert_valid(Plus,'direct positive difference');
    assert_valid(Minus,'direct negative difference');
    DirectFD(iStep).step = h;
    DirectFD(iStep).drTorDx = ([Plus.rTor]-[Minus.rTor])./([Plus.x]-[Minus.x]);
    DirectFD(iStep).drSaarelmaDx = ([Plus.rSaarelma]-[Minus.rSaarelma])./([Plus.x]-[Minus.x]);
    DirectFD(iStep).dPhiTorDpsiN = ([Plus.phiTor]-[Minus.phiTor])./(2*h);
    DirectFD(iStep).dVdx = ([Plus.volume]-[Minus.volume])./([Plus.x]-[Minus.x]);
    DirectFD(iStep).dVdrTor = ([Plus.volume]-[Minus.volume])./([Plus.rTor]-[Minus.rTor]);
    DirectFD(iStep).dVdrSaarelma = ([Plus.volume]-[Minus.volume])./([Plus.rSaarelma]-[Minus.rSaarelma]);
end

% Differentiate the returned lookup coordinates with a small centered step.
% These checks measure interpolation consistency, whereas DirectFD above
% measures the agreement of area integrals with boundary-integral slopes.
h = 1e-7;
Center = lookup_flux_surface_labels(g,queryPsi);
Plus = lookup_flux_surface_labels(g,queryPsi+h);
Minus = lookup_flux_surface_labels(g,queryPsi-h);
assert_valid(Plus,'cached positive difference');
assert_valid(Minus,'cached negative difference');
lookupDerivativeError = [relative_error(([Plus.rTor]-[Minus.rTor])./([Plus.x]-[Minus.x]),[Center.drTorDx]); relative_error(([Plus.rSaarelma]-[Minus.rSaarelma])./([Plus.x]-[Minus.x]),[Center.drSaarelmaDx]); relative_error(([Plus.phiTor]-[Minus.phiTor])/(2*h),[Center.dPhiTorDpsiN])];
assert(max(lookupDerivativeError(:)) < 1e-4,'Lookup derivatives disagree with returned coordinates');

% Each radius y obeys Vy=Vx/(dy/dx), My=Mx*(dy/dx)^2 and
% Ky=Kx*(dy/dx), where K=Vprime*M. Check direct and cached outputs.
transformError = max([transformation_error(Direct),transformation_error(Cached),transformation_error(Fine)]);
assert(transformError < 1e-12,'Coordinate transformation identities failed');
fields = {'x','volume','phiTor','rTor','rSaarelma','dPhiTorDpsiN','drTorDx','drSaarelmaDx','dVdx','dVdrTor','dVdrSaarelma','gradXSqAvg','gradRTorSqAvg','gradRSaarelmaSqAvg','transportMetric','transportMetricRTor','transportMetricRSaarelma'};
Results = struct();
Results.gfileName = gfileName;
Results.table = Table;
Results.psiN = queryPsi;
Results.fixedPsiN = psiN;
Results.Direct = Direct;
Results.Cached = Cached;
Results.Base = Base;
Results.Fine = Fine;
Results.DirectFD = DirectFD;
Results.lookupDerivativeError = lookupDerivativeError;
Results.transformError = transformError;
fprintf('\nMaximum relative discrepancies on 0.85 <= psiN <= 0.999\n');
fprintf('%-25s %14s %14s\n','Field','Lookup/direct','Fine/base');
for iField = 1:numel(fields)
    fieldName = fields{iField};
    Results.lookupError.(fieldName) = relative_error([Cached.(fieldName)],[Direct.(fieldName)]);
    Results.refinementError.(fieldName) = relative_error([Fine.(fieldName)],[Base.(fieldName)]);
    fprintf('%-25s %14.5g %14.5g\n',fieldName,max(Results.lookupError.(fieldName)),max(Results.refinementError.(fieldName)));
end
derivativeFields = fieldnames(DirectFD);
derivativeFields(strcmp(derivativeFields,'step')) = [];
fprintf('\nDirect finite-difference maximum relative discrepancies (h, h/2)\n');
for iField = 1:numel(derivativeFields)
    fieldName = derivativeFields{iField};
    errors = [relative_error(DirectFD(1).(fieldName),[Fine.(fieldName)]);relative_error(DirectFD(2).(fieldName),[Fine.(fieldName)])];
    Results.directDerivativeError.(fieldName) = errors;
    fprintf('%-25s %14.5g %14.5g\n',fieldName,max(errors(1,:)),max(errors(2,:)));
end
fprintf('Maximum cached derivative residual: %.5g\n',max(lookupDerivativeError(:)));
fprintf('Maximum coordinate transformation residual: %.5g\n',transformError);

% Print actual pedestal values and fractional ranges, rather than imposing
% a universal tolerance on the physical constant-Vprime approximation.
Results.Pedestal = table(psiN.',[Fine.dVdx].',[Fine.dVdrTor].',[Fine.dVdrSaarelma].',[Fine.gradXSqAvg].',[Fine.gradRTorSqAvg].',[Fine.gradRSaarelmaSqAvg].','VariableNames',{'psiN','dVdx','dVdrTor','dVdrSaarelma','gradXSqAvg','gradRTorSqAvg','gradRSaarelmaSqAvg'});
disp(Results.Pedestal);
vFields = {'dVdx','dVdrTor','dVdrSaarelma'};
for upper = [0.99,0.999]
    selected = psiN >= 0.95 & psiN <= upper;
    fprintf('Vprime range/mean over 0.95 <= psiN <= %.3f:',upper);
    for iField = 1:numel(vFields)
        values = [Fine(selected).(vFields{iField})];
        fprintf(' %s %.3f%%',vFields{iField},100*(max(values)-min(values))/mean(values));
    end
    fprintf('\n');
end
Results.Figure = [];
if makePlots
    Results.Figure = plot_coordinates(g,Results);
end
end

function radii = radii_at_flux_levels(g,levels)
% Select the first outboard crossing before the axis-height separatrix.
rLine = linspace(g.rmaxis,g.r(end),max(512,2*numel(g.r)));
psiLine = calc_psiN(g,rLine,g.zmaxis*ones(size(rLine)),1);
iSep = find(psiLine(1:end-1) < 1 & psiLine(2:end) >= 1,1,'first');
assert(~isempty(iSep),'No outboard separatrix crossing');
radii = nan(size(levels));
for iLevel = 1:numel(levels)
    iCross = find(psiLine(1:iSep) <= levels(iLevel) & psiLine(2:iSep+1) >= levels(iLevel),1,'first');
    assert(~isempty(iCross),'No outboard crossing for psiN %.8f',levels(iLevel));
    radii(iLevel) = fzero(@(r) calc_psiN(g,r,g.zmaxis,1)-levels(iLevel),rLine(iCross:iCross+1));
end
end

function assert_valid(Geo,description)
% A missing surface is a diagnostic failure, not a point to drop silently.
for iPoint = 1:numel(Geo)
    assert(Geo(iPoint).iErr == 0,'%s failed at psiN %.8f: %s',description,Geo(iPoint).psiN,Geo(iPoint).errorId);
end
end

function error = relative_error(values,reference)
% Diagnostic surfaces avoid coordinate zeros; relative errors are unscaled.
error = abs(values-reference)./max(abs(reference),realmin);
end

function error = transformation_error(Geo)
% Check chain-rule identities for both physical radial coordinates.
error = 0;
suffixes = {'RTor','RSaarelma'};
slopes = {'drTorDx','drSaarelmaDx'};
volumeFields = {'dVdrTor','dVdrSaarelma'};
for iLabel = 1:2
    slope = [Geo.(slopes{iLabel})];
    volume = [Geo.(volumeFields{iLabel})];
    metric = [Geo.(['grad',suffixes{iLabel},'SqAvg'])];
    product = [Geo.(['transportMetric',suffixes{iLabel}])];
    errors = [relative_error(volume.*slope,[Geo.dVdx]),relative_error(metric,[Geo.gradXSqAvg].*slope.^2),relative_error(product,[Geo.transportMetric].*slope),relative_error(product,volume.*metric)];
    error = max(error,max(errors));
end
end

function Figure = plot_coordinates(g,Results)
% A common psiN axis compares physical coordinates on the same surfaces.
% Lines are cached queries; markers are direct finer-contour calculations.
psiN = linspace(0.85,0.999,301);
radii = radii_at_flux_levels(g,psiN);
Geo = lookup_flux_surface_geometry(g,radii,g.zmaxis);
assert_valid(Geo,'plot lookup');
Figure = figure('Color','w','Name','Flux-surface coordinate verification','Position',[100,100,1200,750]);
tiledlayout(2,3);
fieldGroups = {{'dVdx','dVdrTor','dVdrSaarelma'},{'gradXSqAvg','gradRTorSqAvg','gradRSaarelmaSqAvg'},{'transportMetric','transportMetricRTor','transportMetricRSaarelma'}};
yLabels = {'Vprime (m^2)','<|grad label|^2>','Vprime <|grad label|^2> (m^2)'};
colors = lines(3);
for iPanel = 1:3
    nexttile; hold on;
    fields = fieldGroups{iPanel};
    for iLabel = 1:3
        fieldName = fields{iLabel};
        plot(psiN,[Geo.(fieldName)],'Color',colors(iLabel,:),'LineWidth',1.5);
        plot(Results.fixedPsiN,[Results.Fine.(fieldName)],'o','Color',colors(iLabel,:),'HandleVisibility','off');
    end
    grid on; xlabel('\psi_N'); ylabel(yLabels{iPanel});
    if iPanel == 1
        legend('x_{OMP}','r_{tor}','r_{Saarelma}','Location','best');
    end
end
nexttile;
plot(psiN,[Geo.rTor],psiN,[Geo.rSaarelma],'LineWidth',1.5);
grid on; xlabel('\psi_N'); ylabel('Radius (m)'); legend('r_{tor}','r_{Saarelma}','Location','best');
nexttile;
plot(psiN,[Geo.drTorDx],psiN,[Geo.drSaarelmaDx],'LineWidth',1.5);
grid on; xlabel('\psi_N'); ylabel('dr/dx_{OMP}'); legend('r_{tor}','r_{Saarelma}','Location','best');
nexttile; hold on;
for iLabel = 1:3
    fieldName = fieldGroups{1}{iLabel};
    semilogy(Results.psiN,max(Results.lookupError.(fieldName),eps),'Color',colors(iLabel,:),'LineWidth',1.5);
end
set(gca,'YScale','log'); grid on; xlabel('\psi_N'); ylabel('Relative lookup error in Vprime');
sgtitle('Cached geometry (lines) and direct finer contours (markers)');
end
