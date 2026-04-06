clearvars;

% Compare discretized current filaments against the analytic circular-coil
% evaluation for an MPEX configuration.
config_name = 'D1-1';
verbose = 1;
DO_TIMING = 1;

[fil,cur] = setup_MPEX_coils(config_name,verbose);

ntheta_per_wind = 101;  % >= 100 gives reasonable B-component accuracy
ibuild = 0;
for i = 1:fil.ncoils
    if abs(cur(i)) > 1e-8
        [coil0,current0] = build_circular_coil( ...
            fil.rr1(i),fil.rr2(i),fil.z0(i),fil.cl(i), ...
            fil.nturns(i),fil.nlayers(i),cur(i),ntheta_per_wind);
        if ibuild == 0
            coil = coil0;
            current = current0;
        else
            coil = [coil;coil0]; %#ok<AGROW>
            current = [current;current0]; %#ok<AGROW>
        end
        ibuild = ibuild + 1;
    end
end
if verbose
    fprintf('Built %d out of %d coils for the discretized filament test.\n',ibuild,length(cur))
end

bfield.coil = coil;
bfield.current = current;
bfield.type = 'just_coils';

[coil2,current2] = build_MPEX_coils_jackson(config_name,verbose);
bfield2.coil = coil2;
bfield2.current = current2;
bfield2.type = 'MPEX';

nowarn = 0;
npts = 5;
Z = linspace(0,5,npts);
R = linspace(0.001,.15,npts);
phi_radian = linspace(-pi,pi,npts);

% phi_radian = 0;
% R = 0.05;
% Z = 1.05;
filament_eval = @() bfield_general_rzphi(R,Z,phi_radian,bfield,nowarn);
analytic_eval = @() bfield_general_rzphi(R,Z,phi_radian,bfield2,nowarn);

[Bout,ierr] = filament_eval();
[Bout2,ierr2] = analytic_eval();

rel_err = sqrt((Bout.br-Bout2.br).^2 + (Bout.bphi-Bout2.bphi).^2 + (Bout.bz-Bout2.bz).^2) ...
    ./ sqrt(Bout2.br.^2 + Bout2.bz.^2) * 100;

fprintf('\n')
fprintf('===============================================================\n')
fprintf('MPEX discretization test for configuration %s\n',config_name)
fprintf('---------------------------------------------------------------\n')
fprintf('Evaluation points : %d\n',length(R))
fprintf('Max rel. error    : %.6f %%\n',max(rel_err))
fprintf('Mean rel. error   : %.6f %%\n',mean(rel_err))
fprintf('ierr flags        : filaments=%d analytic=%d\n',ierr,ierr2)
if DO_TIMING
    t_filaments = timeit(filament_eval);
    t_analytic = timeit(analytic_eval);
    fprintf('Filament timing   : %.6f s (timeit)\n',t_filaments)
    fprintf('Analytic timing   : %.6f s (timeit)\n',t_analytic)
    fprintf('Timing ratio      : %.2fx\n',t_filaments/t_analytic)
end
fprintf('===============================================================\n')

if length(R) <= 10
    fprintf('\n')
    fprintf('Pointwise comparison\n')
    fprintf('-------------------------------------------------------------------------------\n')
    for i = 1:length(R)
        fprintf('Point %d at [R,Z,phi] = [%8.5f, %8.5f, %8.5f]\n',i,R(i),Z(i),phi_radian(i))
        fprintf('  FILAMENTS : [Br,Bz,Bphi] = [% .6e, % .6e, % .6e]\n',Bout.br(i),Bout.bz(i),Bout.bphi(i))
        fprintf('  ANALYTIC  : [Br,Bz,Bphi] = [% .6e, % .6e, % .6e]\n',Bout2.br(i),Bout2.bz(i),Bout2.bphi(i))
        fprintf('  Rel. error: %.6f %%\n',rel_err(i))
    end
end
