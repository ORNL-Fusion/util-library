clearvars;

gfile_name = '/Users/jjl/ORNL Dropbox/Jeremy Lore/DIII-D/Haskey/g179447.02280';
g = readg_g3d(gfile_name);

compile_calc_psi_mex();

TIMEIT = 1;
PLOTIT = 0;
RNG_SEED = 1;
NPTS_COMPARE = 10000;

rng(RNG_SEED);

Rmin = g.rgrid1;
Rmax = g.rgrid1 + g.xdim;
Zmin = g.zmid - 0.5*g.zdim;
Zmax = g.zmid + 0.5*g.zdim;

Rcmp = Rmin + (Rmax - Rmin)*rand(1,NPTS_COMPARE);
Zcmp = Zmin + (Zmax - Zmin)*rand(1,NPTS_COMPARE);

[psi0,ierr0,dpsidr0,dpsidz0] = calc_psi(g,Rcmp,Zcmp,1);
[psi1,ierr1,dpsidr1,dpsidz1] = calc_psi_mex(g,Rcmp,Zcmp,1);
[B0,ierrB0] = bfield_geq_bicub(g,Rcmp,Zcmp,1);
[B1,ierrB1] = bfield_geq_bicub_mex(g,Rcmp,Zcmp,1);

assert(isequal(ierr0,ierr1),'calc_psi ierr mismatch')
assert(isequal(ierrB0,ierrB1),'bfield_geq_bicub ierr mismatch')
assert(max(abs(psi0(:) - psi1(:)),[],'omitnan') < 2e-12,'psi mismatch')
assert(max(abs(dpsidr0(:) - dpsidr1(:)),[],'omitnan') < 2e-12,'dpsidr mismatch')
assert(max(abs(dpsidz0(:) - dpsidz1(:)),[],'omitnan') < 2e-12,'dpsidz mismatch')
assert(max(abs(B0.br(:) - B1.br(:)),[],'omitnan') < 2e-12,'Br mismatch')
assert(max(abs(B0.bz(:) - B1.bz(:)),[],'omitnan') < 2e-12,'Bz mismatch')
assert(max(abs(B0.bphi(:) - B1.bphi(:)),[],'omitnan') < 2e-12,'Bphi mismatch')

fprintf('\n')
fprintf('bfield_geq_bicub mex timing study\n')
fprintf('gfile: %s\n',gfile_name)
fprintf('comparison points      : %d\n',NPTS_COMPARE)

bench_many_points.name = 'many point bulk eval';
bench_many_points.R = Rmin + (Rmax - Rmin)*rand(1,1000);
bench_many_points.Z = Zmin + (Zmax - Zmin)*rand(1,1000);
bench_many_points.nreps = 1;

bench_many_steps.name = 'small vector repeated eval';
bench_many_steps.R = Rmin + (Rmax - Rmin)*rand(1,4);
bench_many_steps.Z = Zmin + (Zmax - Zmin)*rand(1,4);
bench_many_steps.nreps = 4000;

if TIMEIT
    compare_timing_case(g,bench_many_points);
    compare_timing_case(g,bench_many_steps);
end

if PLOTIT
    figure
    subplot(1,2,1)
    histogram(B0.bphi - B1.bphi,100)
    title('B_\phi difference')
    subplot(1,2,2)
    histogram(psi0 - psi1,100)
    title('\psi difference')
end

function compare_timing_case(g,bench)
f_orig = @() run_bfield_case(g,bench,@bfield_geq_bicub);
f_mex = @() run_bfield_case(g,bench,@bfield_geq_bicub_mex);

elapsed_orig = time_this(f_orig);
elapsed_mex = time_this(f_mex);

nevals = numel(bench.R) * bench.nreps;

fprintf('\n')
fprintf('Benchmark: %s\n',bench.name)
fprintf('  nPoints              : %d\n',numel(bench.R))
fprintf('  repetitions          : %d\n',bench.nreps)
fprintf('  original elapsed [s] : %.6e\n',elapsed_orig)
fprintf('  original evals/s     : %.6e\n',nevals/elapsed_orig)
fprintf('  mex elapsed [s]      : %.6e\n',elapsed_mex)
fprintf('  mex evals/s          : %.6e\n',nevals/elapsed_mex)
fprintf('  speedup old/mex      : %.3fx\n',elapsed_orig/elapsed_mex)
end

function run_bfield_case(g,bench,bfield_fun)
for i = 1:bench.nreps
    [~,ierr] = bfield_fun(g,bench.R,bench.Z,1);
    assert(all(~ierr),'Encountered off-grid point during timing case')
end
end

function elapsed = time_this(fun_handle)
fun_handle();

if exist('timeit','file') == 2
    elapsed = timeit(fun_handle);
else
    nruns = 5;
    tic;
    for k = 1:nruns
        fun_handle();
    end
    elapsed = toc / nruns;
end
end
