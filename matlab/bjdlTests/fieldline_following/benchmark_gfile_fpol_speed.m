function Results = benchmark_gfile_fpol_speed(gfileName)
% Compare cached polynomial and PCHIP F evaluation on a real gfile.
% Results = benchmark_gfile_fpol_speed(gfileName)
% The default gfile is the MAST-U case under C:\Work. Timings exclude
% polynomial fitting and PCHIP construction, which happen once per gfile.

if nargin < 1 || isempty(gfileName)
    gfileName = 'C:\Work\SOLPS\runs\mastu\cd-51685-400ms\baserun\gfile';
end
testDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(testDirectory,'..','..','bfield_library_jdl','EFIT_routines'));
addpath(fullfile(testDirectory,'..','..','bfield_library_jdl','calc_B'));
addpath(fullfile(testDirectory,'..','..','bfield_library_jdl','util'));
g = readg_g3d(gfileName);
if isempty(g)
    error('benchmark_gfile_fpol_speed:NoGfile','Could not read %s',gfileName)
end
polyCoeffs = polyfit(g.pn,g.fpol,7);

% Select reproducible core points across the outboard axis-height line.
% All timed calls use the same points and the same precomputed coefficients.
rLine = linspace(g.rmaxis,g.r(end),max(512,2*numel(g.r)));
psiNLine = calc_psiN(g,rLine,g.zmaxis*ones(size(rLine)),1);
iSep = find(psiNLine(1:end-1) < 1 & psiNLine(2:end) >= 1,1,'first');
if isempty(iSep)
    error('benchmark_gfile_fpol_speed:NoSeparatrix','No outboard separatrix crossing was found')
end
rSep = fzero(@(rValue) calc_psiN(g,rValue,g.zmaxis,1)-1,rLine(iSep:iSep+1));
caseSizes = [1,4,32,1000,10000];
rngState = rng;
rng(1);
rPool = g.rmaxis+(rSep-g.rmaxis)*(0.02+0.96*rand(1,2*max(caseSizes)));
rng(rngState);
psiNPool = calc_psiN(g,rPool,g.zmaxis*ones(size(rPool)),1);
rPool = rPool(isfinite(psiNPool) & psiNPool > 0 & psiNPool < 1);
if numel(rPool) < max(caseSizes)
    error('benchmark_gfile_fpol_speed:TooFewPoints','Could not find enough core points for the benchmark')
end

% Check that the local PCHIP path reproduces the production field routine.
% The polynomial differs in Bphi by its fitting error; Br and Bz are common.
rCheck = rPool(1:1000);
zCheck = g.zmaxis*ones(size(rCheck));
[Production,ierr] = bfield_geq_bicub(g,rCheck,zCheck,1);
PchipField = evaluate_field(g,rCheck,zCheck,polyCoeffs,true);
PolyField = evaluate_field(g,rCheck,zCheck,polyCoeffs,false);
assert(~any(ierr) && max(abs(Production.br-PchipField.br)) < 1e-12);
assert(max(abs(Production.bz-PchipField.bz)) < 1e-12);
assert(max(abs(Production.bphi-PchipField.bphi)) < 1e-12);
fprintf('\nF evaluation benchmark: %s\n',gfileName)
fprintf('Polynomial node error: %.3g T m\n',max(abs(polyval(polyCoeffs,g.pn)-g.fpol)))
fprintf('Maximum Bphi difference on 1000 core points: %.3g T\n',max(abs(PolyField.bphi-PchipField.bphi)))

% Time F alone and the full field calculation at common call sizes.
% Warm up each path, alternate timing order, and report the median of five
% trials. Ratios greater than one mean PCHIP is slower than the polynomial.
Results = struct();
Results.gfileName = gfileName;
Results.nPoints = caseSizes;
Results.repetitions = [10000,10000,5000,500,50];
Results.trialSeconds = nan(4,5,numel(caseSizes));
Results.fpolPolySeconds = nan(size(caseSizes));
Results.fpolPchipSeconds = nan(size(caseSizes));
Results.fieldPolySeconds = nan(size(caseSizes));
Results.fieldPchipSeconds = nan(size(caseSizes));
fprintf('\n%8s  %12s  %12s  %9s  %12s  %12s  %9s\n','points','F poly [us]','F PCHIP [us]','ratio','field poly [us]','field PCHIP [us]','ratio')
for iCase = 1:numel(caseSizes)
    nPoints = caseSizes(iCase);
    rPoints = rPool(1:nPoints);
    zPoints = g.zmaxis*ones(size(rPoints));
    psiNPoints = calc_psiN(g,rPoints,zPoints,1);
    methods = cell(1,4);
    methods{1} = @() evaluate_fpol(g,psiNPoints,polyCoeffs,false);
    methods{2} = @() evaluate_fpol(g,psiNPoints,polyCoeffs,true);
    methods{3} = @() evaluate_field(g,rPoints,zPoints,polyCoeffs,false);
    methods{4} = @() evaluate_field(g,rPoints,zPoints,polyCoeffs,true);
    for iMethod = 1:4
        for iWarmup = 1:20
            methods{iMethod}();
        end
    end
    for iTrial = 1:5
        if mod(iTrial,2) == 1
            order = [1,2,3,4];
        else
            order = [2,1,4,3];
        end
        for iMethod = order
            evaluate = methods{iMethod};
            tic;
            for iRepeat = 1:Results.repetitions(iCase)
                lastValue = evaluate();
            end
            Results.trialSeconds(iMethod,iTrial,iCase) = toc/Results.repetitions(iCase);
            if isempty(lastValue)
                error('benchmark_gfile_fpol_speed:EmptyResult','A benchmark evaluator returned no data')
            end
        end
    end
    medianSeconds = median(Results.trialSeconds(:,:,iCase),2);
    Results.fpolPolySeconds(iCase) = medianSeconds(1);
    Results.fpolPchipSeconds(iCase) = medianSeconds(2);
    Results.fieldPolySeconds(iCase) = medianSeconds(3);
    Results.fieldPchipSeconds(iCase) = medianSeconds(4);
    fpolPolyUs = 1e6*Results.fpolPolySeconds(iCase);
    fpolPchipUs = 1e6*Results.fpolPchipSeconds(iCase);
    fieldPolyUs = 1e6*Results.fieldPolySeconds(iCase);
    fieldPchipUs = 1e6*Results.fieldPchipSeconds(iCase);
    fpolRatio = fpolPchipUs/fpolPolyUs;
    fieldRatio = fieldPchipUs/fieldPolyUs;
    fprintf('%8d  %12.3f  %12.3f  %9.3f  %15.3f  %16.3f  %9.3f\n',nPoints,fpolPolyUs,fpolPchipUs,fpolRatio,fieldPolyUs,fieldPchipUs,fieldRatio)
end
end

function fpol = evaluate_fpol(g,psiN,polyCoeffs,usePchip)
% Change only the one-dimensional F evaluator between benchmark cases.
if usePchip
    fpol = ppval(g.fpol_pp,psiN);
else
    fpol = polyval(polyCoeffs,psiN);
end
end

function Bout = evaluate_field(g,R,Z,polyCoeffs,usePchip)
% Match the core branch of bfield_geq_bicub with either F evaluator.
[psi,ierr,dpsidr,dpsidz] = calc_psi(g,R,Z,1);
rInv = 1./R;
psiN = (g.ip_sign*psi-g.ssimag)/(g.ssibry-g.ssimag);
Bout.br = -dpsidz.*rInv;
Bout.bz = dpsidr.*rInv;
Bout.bphi = evaluate_fpol(g,psiN,polyCoeffs,usePchip).*rInv;
Bout.br(ierr) = nan;
Bout.bz(ierr) = nan;
Bout.bphi(ierr) = nan;
end
