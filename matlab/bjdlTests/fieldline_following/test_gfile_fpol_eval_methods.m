clearvars;

gfile_name = '/Users/jjl/ORNL Dropbox/Jeremy Lore/DIII-D/Haskey/g179447.02280';
g = readg_g3d(gfile_name);

PLOTIT = 0;
TIMEIT = 1;
NPTS = 100000;
RNG_SEED = 1;

rng(RNG_SEED);

Rmin = g.rgrid1;
Rmax = g.rgrid1 + g.xdim;
Zmin = g.zmid - 0.5*g.zdim;
Zmax = g.zmid + 0.5*g.zdim;

R = Rmin + (Rmax - Rmin)*rand(1,NPTS);
Z = Zmin + (Zmax - Zmin)*rand(1,NPTS);

[psiN,~,ierr] = calc_psiN(g,R,Z,1);

is_valid = ~ierr & isfinite(psiN);
is_inside = is_valid & psiN <= 1;

psiN_inside = psiN(is_inside);
R_inside = R(is_inside);

fpol_poly_inside = polyval(g.fpol_coeffs,psiN_inside);

fpol_linear_inside = interp1(g.pn,g.fpol,psiN_inside,'linear','extrap');

fpol_pchip_inside = interp1(g.pn,g.fpol,psiN_inside,'pchip','extrap');

fpol_uniform_inside = interp_uniform_linear(g.pn,g.fpol,psiN_inside);

bphi_poly_inside = fpol_poly_inside./R_inside;
bphi_linear_inside = fpol_linear_inside./R_inside;
bphi_pchip_inside = fpol_pchip_inside./R_inside;
bphi_uniform_inside = fpol_uniform_inside./R_inside;

fprintf('\n')
fprintf('gfile fpol evaluation comparison\n')
fprintf('gfile: %s\n',gfile_name)
fprintf('random points sampled : %d\n',NPTS)
fprintf('valid in-grid points  : %d\n',nnz(is_valid))
fprintf('inside separatrix     : %d\n',nnz(is_inside))

print_diff_stats('polyfit(7) vs linear interp: fpol inside sep', ...
    fpol_poly_inside,fpol_linear_inside)
print_diff_stats('polyfit(7) vs pchip interp: fpol inside sep', ...
    fpol_poly_inside,fpol_pchip_inside)
print_diff_stats('linear interp vs uniform-linear: fpol inside sep', ...
    fpol_linear_inside,fpol_uniform_inside)

print_diff_stats('polyfit(7) vs linear interp: Bphi inside sep', ...
    bphi_poly_inside,bphi_linear_inside)
print_diff_stats('polyfit(7) vs pchip interp: Bphi inside sep', ...
    bphi_poly_inside,bphi_pchip_inside)
print_diff_stats('linear interp vs uniform-linear: Bphi inside sep', ...
    bphi_linear_inside,bphi_uniform_inside)

if TIMEIT
    fprintf('\n')
    fprintf('Timing on inside-separatrix psiN sample\n')

    time_eval_method('polyfit(7) / polyval', ...
        @() polyval(g.fpol_coeffs,psiN_inside),numel(psiN_inside))
    time_eval_method('interp1 linear', ...
        @() interp1(g.pn,g.fpol,psiN_inside,'linear','extrap'),numel(psiN_inside))
    time_eval_method('interp1 pchip', ...
        @() interp1(g.pn,g.fpol,psiN_inside,'pchip','extrap'),numel(psiN_inside))
    time_eval_method('uniform linear', ...
        @() interp_uniform_linear(g.pn,g.fpol,psiN_inside),numel(psiN_inside))
end

if PLOTIT
    figure
    tiledlayout(2,2)

    nexttile
    plot(psiN_inside,fpol_poly_inside - fpol_linear_inside,'.','markersize',4)
    xlabel('\psi_N')
    ylabel('\Delta F [T m]')
    title('polyfit(7) - linear interp')

    nexttile
    plot(psiN_inside,fpol_poly_inside - fpol_pchip_inside,'.','markersize',4)
    xlabel('\psi_N')
    ylabel('\Delta F [T m]')
    title('polyfit(7) - pchip interp')

    nexttile
    histogram(fpol_poly_inside - fpol_linear_inside,100)
    xlabel('\Delta F [T m]')
    ylabel('count')
    title('polyfit(7) - linear interp')

    nexttile
    histogram(bphi_poly_inside - bphi_linear_inside,100)
    xlabel('\Delta B_\phi [T]')
    ylabel('count')
    title('B_\phi: polyfit(7) - linear interp')
end

function print_diff_stats(label,a,b)
da = a - b;
abs_ref = max(abs(b),eps);

fprintf('\n')
fprintf('%s\n',label)
fprintf('  mean abs diff        : %.6e\n',mean(abs(da)))
fprintf('  rms diff             : %.6e\n',sqrt(mean(da.^2)))
fprintf('  max abs diff         : %.6e\n',max(abs(da)))
fprintf('  mean rel diff        : %.6e\n',mean(abs(da)./abs_ref))
fprintf('  max rel diff         : %.6e\n',max(abs(da)./abs_ref))
end

function time_eval_method(label,fun_handle,npts)
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

fprintf('%s\n',label)
fprintf('  elapsed/call [s]     : %.6e\n',elapsed)
fprintf('  evaluations/s        : %.6e\n',npts/elapsed)
end

function yq = interp_uniform_linear(x,y,xq)
nx = length(x);
x0 = x(1);
dx = x(2) - x(1);

t = (xq - x0)./dx + 1;
i0 = floor(t);
w = t - i0;

left = i0 < 1;
right = i0 >= nx;

i0(left) = 1;
i0(right) = nx - 1;
w(left) = t(left) - 1;
w(right) = t(right) - (nx - 1);

i0 = max(1,min(nx-1,i0));
i1 = i0 + 1;

yq = (1 - w).*y(i0) + w.*y(i1);
end
