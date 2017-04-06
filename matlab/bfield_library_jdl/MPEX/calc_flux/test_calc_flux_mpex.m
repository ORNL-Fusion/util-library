clearvars;


% Set up magnetic field -- MPEX parameters
helicon_current = 440;
current_A = 6400;
current_B = 6400;
current_C = [];
config = 'flat';
verbose = 0;

% Build bfield structure
[coil,current] = build_Proto_coils_jackson(helicon_current,current_A,current_B,config,verbose,current_C);
bfield.coil = coil; bfield.current = current; bfield.type = 'MPEX';

% Evaluate flux at target

% % First define geometry
skimmer = 1;
plasma_radius_cm = 1;
target_position = 2;
sleeve = 1;
plot_geo = 1; if plot_geo; newfig_geo = 1; else; newfig_geo = 0; end
geo = get_Proto_geometry(plot_geo,newfig_geo,skimmer,target_position,sleeve);

% Interpolate along target radius
RMAX = 0.1221; % Inner radius of coil
ZTARG = 4.33415; % Axial position of target
ninterp = 1000;
zinterp = ZTARG*ones(1,ninterp);
rinterp = linspace(0,RMAX,ninterp);
psiinterp = calc_psi_mpex(coil,current,rinterp,zinterp);

% figure; hold on; box on;
% plot(rinterp,psiinterp)

% Evaluation and mapping to target
reval = 0.05;
zeval = 1.05;
psi_eval = calc_psi_mpex(coil,current,reval,zeval);
rmap = interp1(psiinterp,rinterp,psi_eval)
