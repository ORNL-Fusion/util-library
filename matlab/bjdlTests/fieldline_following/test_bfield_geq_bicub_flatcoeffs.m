clearvars;

gfile_name = '/Users/jjl/ORNL Dropbox/Jeremy Lore/DIII-D/Haskey/g179447.02280';
g = readg_g3d(gfile_name);
g_flat = add_psi_bicub_coeffs_inv_flat(g);

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
[psi1,ierr1,dpsidr1,dpsidz1] = calc_psi_flatcoeffs(g_flat,Rcmp,Zcmp,1);
[psi2,ierr2,dpsidr2,dpsidz2] = calc_psi_flatcoeffs_horner(g_flat,Rcmp,Zcmp,1);
[B0,ierrB0] = bfield_geq_bicub(g,Rcmp,Zcmp,1);
[B1,ierrB1] = bfield_geq_bicub_flatcoeffs(g_flat,Rcmp,Zcmp,1);
[B2,ierrB2] = bfield_geq_bicub_flatcoeffs_horner(g_flat,Rcmp,Zcmp,1);

assert(isequal(ierr0,ierr1),'calc_psi ierr mismatch')
assert(isequal(ierr0,ierr2),'calc_psi horner ierr mismatch')
assert(isequal(ierrB0,ierrB1),'bfield_geq_bicub ierr mismatch')
assert(isequal(ierrB0,ierrB2),'bfield_geq_bicub horner ierr mismatch')
assert(max(abs(psi0(:) - psi1(:)),[],'omitnan') < 1e-12,'psi mismatch')
assert(max(abs(psi0(:) - psi2(:)),[],'omitnan') < 1e-12,'psi horner mismatch')
assert(max(abs(dpsidr0(:) - dpsidr1(:)),[],'omitnan') < 1e-12,'dpsidr mismatch')
assert(max(abs(dpsidr0(:) - dpsidr2(:)),[],'omitnan') < 1e-12,'dpsidr horner mismatch')
assert(max(abs(dpsidz0(:) - dpsidz1(:)),[],'omitnan') < 1e-12,'dpsidz mismatch')
assert(max(abs(dpsidz0(:) - dpsidz2(:)),[],'omitnan') < 1e-12,'dpsidz horner mismatch')
assert(max(abs(B0.br(:) - B1.br(:)),[],'omitnan') < 1e-12,'Br mismatch')
assert(max(abs(B0.br(:) - B2.br(:)),[],'omitnan') < 1e-12,'Br horner mismatch')
assert(max(abs(B0.bz(:) - B1.bz(:)),[],'omitnan') < 1e-12,'Bz mismatch')
assert(max(abs(B0.bz(:) - B2.bz(:)),[],'omitnan') < 1e-12,'Bz horner mismatch')
assert(max(abs(B0.bphi(:) - B1.bphi(:)),[],'omitnan') < 1e-12,'Bphi mismatch')
assert(max(abs(B0.bphi(:) - B2.bphi(:)),[],'omitnan') < 1e-12,'Bphi horner mismatch')

fprintf('\n')
fprintf('bfield_geq_bicub flat-coeff timing study\n')
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
    compare_timing_case(g,g_flat,bench_many_points);
    compare_timing_case(g,g_flat,bench_many_steps);
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

function compare_timing_case(g,g_flat,bench)
f_orig = @() run_bfield_case(g,bench,@bfield_geq_bicub);
f_flat = @() run_bfield_case(g_flat,bench,@bfield_geq_bicub_flatcoeffs);
f_horner = @() run_bfield_case(g_flat,bench,@bfield_geq_bicub_flatcoeffs_horner);

elapsed_orig = time_this(f_orig);
elapsed_flat = time_this(f_flat);
elapsed_horner = time_this(f_horner);

nevals = numel(bench.R) * bench.nreps;

fprintf('\n')
fprintf('Benchmark: %s\n',bench.name)
fprintf('  nPoints              : %d\n',numel(bench.R))
fprintf('  repetitions          : %d\n',bench.nreps)
fprintf('  original elapsed [s] : %.6e\n',elapsed_orig)
fprintf('  original evals/s     : %.6e\n',nevals/elapsed_orig)
fprintf('  flat elapsed [s]     : %.6e\n',elapsed_flat)
fprintf('  flat evals/s         : %.6e\n',nevals/elapsed_flat)
fprintf('  speedup old/flat     : %.3fx\n',elapsed_orig/elapsed_flat)
fprintf('  horner elapsed [s]   : %.6e\n',elapsed_horner)
fprintf('  horner evals/s       : %.6e\n',nevals/elapsed_horner)
fprintf('  speedup old/horner   : %.3fx\n',elapsed_orig/elapsed_horner)
fprintf('  speedup flat/horner  : %.3fx\n',elapsed_flat/elapsed_horner)
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
