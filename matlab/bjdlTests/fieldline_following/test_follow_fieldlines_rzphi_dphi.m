clearvars;

gfile_name = '/Users/jjl/ORNL Dropbox/Jeremy Lore/DIII-D/Haskey/g179447.02280';
g = readg_g3d(gfile_name);

PLOTIT = 0;
TIMEIT = 1;
TEST = 'COMPARE';  % 'A' = original, 'B' = direct dphi mex worker, 'COMPARE' = both

Bfield.g = g;
Bfield.type = 'gfile';
Bfield.nsym = 1;

%% Ensure mex exists before testing mex path
compile_calc_psi_mex();

%% Baseline correctness test
Rstart = 1.8;
Rend = 2.0;
nTest = 10;
RTest = linspace(Rstart,Rend,nTest);
ZTest = zeros(size(RTest));
phiTest = 0;
nsteps = 4000;
dphi = -0.1;

[f,ierr,i_last_good] = follow_fieldlines_rzphi_dphi(Bfield,RTest,ZTest,phiTest,dphi,nsteps);
[f_mex,ierr_mex,i_last_good_mex] = follow_fieldlines_dphi_gfile_mex(Bfield,RTest,ZTest,phiTest,dphi,nsteps);

if PLOTIT
    plot_gfile(g)
    plot(f.r,f.z,'m','LineWidth',3)
end

assert(ierr == 0,'Fieldline following returned ierr = %d',ierr)
assert(i_last_good == nsteps + 1,'Unexpected i_last_good = %d',i_last_good)
assert(ierr_mex == 0,'MEX fieldline following returned ierr = %d',ierr_mex)
assert(i_last_good_mex == nsteps + 1,'Unexpected mex i_last_good = %d',i_last_good_mex)
assert(isequal(size(f.r),[nsteps+1,nTest]),'Unexpected size for f.r')
assert(isequal(size(f.z),[nsteps+1,nTest]),'Unexpected size for f.z')
assert(isequal(size(f.phi),[nsteps+1,1]),'Unexpected size for f.phi')
assert(isequal(size(f_mex.r),[nsteps+1,nTest]),'Unexpected size for f_mex.r')
assert(isequal(size(f_mex.z),[nsteps+1,nTest]),'Unexpected size for f_mex.z')
assert(isequal(size(f_mex.phi),[nsteps+1,1]),'Unexpected size for f_mex.phi')
assert(max(abs(f.r(1,:) - RTest)) < 1e-12,'Initial R values were not preserved')
assert(max(abs(f.z(1,:) - ZTest)) < 1e-12,'Initial Z values were not preserved')
assert(abs(f.phi(1) - phiTest) < 1e-12,'Initial phi value was not preserved')
assert(max(abs(f_mex.r(1,:) - RTest)) < 1e-12,'Initial mex R values were not preserved')
assert(max(abs(f_mex.z(1,:) - ZTest)) < 1e-12,'Initial mex Z values were not preserved')
assert(abs(f_mex.phi(1) - phiTest) < 1e-12,'Initial mex phi value was not preserved')
assert(all(isfinite(f.r(:))),'Encountered non-finite values in f.r')
assert(all(isfinite(f.z(:))),'Encountered non-finite values in f.z')
assert(all(isfinite(f.phi(:))),'Encountered non-finite values in f.phi')
assert(all(isfinite(f_mex.r(:))),'Encountered non-finite values in f_mex.r')
assert(all(isfinite(f_mex.z(:))),'Encountered non-finite values in f_mex.z')
assert(all(isfinite(f_mex.phi(:))),'Encountered non-finite values in f_mex.phi')
assert(max(abs(f.r(:) - f_mex.r(:))) < 2e-12,'MEX and generic f.r do not match')
assert(max(abs(f.z(:) - f_mex.z(:))) < 2e-12,'MEX and generic f.z do not match')
assert(max(abs(f.phi(:) - f_mex.phi(:))) < 1e-12,'MEX and generic f.phi do not match')

%% Timing / benchmarking
benchmarks = build_benchmarks();

if TIMEIT || strcmp(TEST,'COMPARE')
    fprintf('\n')
    fprintf('follow_fieldlines_rzphi_dphi timing study\n')
    fprintf('gfile: %s\n',gfile_name)
end

for ibench = 1:numel(benchmarks)
    bench = benchmarks(ibench);
    if isempty(bench.RTest)
        fprintf('\n')
        fprintf('Skipping "%s" benchmark.\n',bench.name)
        fprintf('Provide closed-surface RTest values here to enable long-step timing.\n')
        continue;
    end

    switch upper(TEST)
        case 'A'
            run_selected_test(Bfield,bench,@follow_fieldlines_rzphi_dphi,'generic dispatcher',TIMEIT)
        case 'B'
            run_selected_test(Bfield,bench,@follow_fieldlines_dphi_gfile_mex,'direct dphi mex worker',TIMEIT)
        case 'COMPARE'
            compare_timing_case(Bfield,bench,TIMEIT)
        otherwise
            error('TEST must be ''A'', ''B'', or ''COMPARE''')
    end
end

function benchmarks = build_benchmarks()
benchmarks = struct([]);

benchmarks(1).name = 'many start points';
benchmarks(1).RTest = linspace(1.8,2.28,1000);
benchmarks(1).ZTest = zeros(size(benchmarks(1).RTest));
benchmarks(1).phiTest = 0;
benchmarks(1).dphi = -0.05;
benchmarks(1).nsteps = 100;
benchmarks(1).require_success = 1;

benchmarks(2).name = 'many integration steps';
benchmarks(2).RTest = [1.7, 1.8, 1.9, 2.0];
benchmarks(2).ZTest = zeros(size(benchmarks(2).RTest));
benchmarks(2).phiTest = 0;
benchmarks(2).dphi = -0.01;
benchmarks(2).nsteps = 4000;
benchmarks(2).require_success = 1;
end

function run_selected_test(Bfield,bench,follow_fun,variant_name,do_timeit)
[~,ierr,i_last_good] = run_case(Bfield,bench,follow_fun);
fprintf('\n')
fprintf('Benchmark: %s\n',bench.name)
fprintf('  variant              : %s\n',variant_name)
fprintf('  ierr                 : %d\n',ierr)
fprintf('  i_last_good          : %d\n',i_last_good)
fprintf('  nStartPoints         : %d\n',numel(bench.RTest))
fprintf('  nsteps               : %d\n',bench.nsteps)
fprintf('  dphi [rad]           : %.6f\n',bench.dphi)

if do_timeit
    elapsed = time_follow_case(Bfield,bench,follow_fun);
    nfieldline_steps = numel(bench.RTest) * bench.nsteps;
    fprintf('  elapsed/call [s]     : %.6e\n',elapsed)
    fprintf('  fieldline-steps/s    : %.6e\n',nfieldline_steps/elapsed)
end
end

function compare_timing_case(Bfield,bench,do_timeit)
[~,ierr_generic,i_last_good_generic] = run_case(Bfield,bench,@follow_fieldlines_rzphi_dphi);
[~,ierr_mex,i_last_good_mex] = run_case(Bfield,bench,@follow_fieldlines_dphi_gfile_mex);

fprintf('\n')
fprintf('Benchmark: %s\n',bench.name)
fprintf('  generic ierr         : %d\n',ierr_generic)
fprintf('  generic i_last_good  : %d\n',i_last_good_generic)
fprintf('  mex ierr             : %d\n',ierr_mex)
fprintf('  mex i_last_good      : %d\n',i_last_good_mex)
fprintf('  nStartPoints         : %d\n',numel(bench.RTest))
fprintf('  nsteps               : %d\n',bench.nsteps)
fprintf('  dphi [rad]           : %.6f\n',bench.dphi)

if do_timeit
    elapsed_generic = time_follow_case(Bfield,bench,@follow_fieldlines_rzphi_dphi);
    elapsed_mex = time_follow_case(Bfield,bench,@follow_fieldlines_dphi_gfile_mex);

    nfieldline_steps = numel(bench.RTest) * bench.nsteps;
    rate_generic = nfieldline_steps / elapsed_generic;
    rate_mex = nfieldline_steps / elapsed_mex;
    speedup_mex = elapsed_generic / elapsed_mex;

    fprintf('  generic elapsed [s]  : %.6e\n',elapsed_generic)
    fprintf('  generic steps/s      : %.6e\n',rate_generic)
    fprintf('  mex elapsed [s]      : %.6e\n',elapsed_mex)
    fprintf('  mex steps/s          : %.6e\n',rate_mex)
    fprintf('  speedup old/mex      : %.3fx\n',speedup_mex)
end
end

function elapsed = time_follow_case(Bfield,bench,follow_fun)
run_once = @() run_case(Bfield,bench,follow_fun);
run_once();

if exist('timeit','file') == 2
    elapsed = timeit(run_once);
else
    nruns = 5;
    tic;
    for k = 1:nruns
        run_once();
    end
    elapsed = toc / nruns;
end
end

function [s,ierr,i_last_good] = run_case(Bfield,bench,follow_fun)
[s,ierr,i_last_good] = follow_fun( ...
    Bfield,bench.RTest,bench.ZTest,bench.phiTest,bench.dphi,bench.nsteps);

if bench.require_success
    assert(ierr == 0, ...
        'Timing case "%s" returned ierr = %d (i_last_good = %d)', ...
        bench.name,ierr,i_last_good)
end
end
