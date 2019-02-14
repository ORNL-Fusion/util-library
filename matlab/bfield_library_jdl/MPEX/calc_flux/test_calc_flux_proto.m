clearvars;

%-----------------------------------------------------
% Set up magnetic field -- MPEX parameters
%-----------------------------------------------------
% Two methods. 
if 1
    % 1) Specify currents corresponding to power supplys and pre-defined configurations
    helicon_current = 200;
    current_A = 6400;
    current_B = 6400;
    current_C = [];
    config = 'flat';
    current_in = [helicon_current,current_A,current_B,current_C];
else
    % 2) Specify a 12 element array corresponding to each coil winding current
    current_in = [6400 0 260 260]; current_in(5:12) = 1200; %current_in(6:7) = 800;
    config = [];
end

% Build coils
verbose = 1;
[coil,current] = build_Proto_coils_jackson(current_in,config,verbose);

%--------------------------------------------------------------------------
% Example 1: 
% Calculate psi at a point R,Z
%--------------------------------------------------------------------------
% Reval = 0.063; Zeval = 1.75;
Reval = 0.045; Zeval = 1.75;
psi_eval = calc_psi_mpex(coil,current,Reval,Zeval);
fprintf('\nExample 1\n')
fprintf('At point (R,Z) = (%6.3f,%6.3f) [m], psi = %6.3e [Wb]\n',Reval,Zeval,psi_eval)

%--------------------------------------------------------------------------
% Example 2: 
% Map a point R,Z to the a line at constant axial position (e.g., target)
% by interpolating axial flux.
%--------------------------------------------------------------------------
Ztarg = 4.33415;  % Set axial position for mapping
% Ztarg = 2.68;  % Set axial position for mapping
Rmap = map_pt_to_target_mpex(coil,current,Ztarg,Reval,Zeval);
fprintf('\nExample 2\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,Rmap,Ztarg)

%--------------------------------------------------------------------------
% Example 3: 
% Map a point R,Z to the a line at constant axial position (e.g., target)
% by following a field line
%--------------------------------------------------------------------------
MAKE_PLOT = 1;

bfield.coil = coil; bfield.current = current; bfield.type = 'MPEX';
dz_want = 0.01;
L = Ztarg - Zeval;
nsteps = round(abs(L/dz_want));
dz = L/nsteps;
f = follow_fieldlines_rzphi_dz(bfield,Reval,Zeval,0,dz,nsteps);
fprintf('\nExample 3\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,f.r(end),Ztarg)

if MAKE_PLOT
    skimmer = 1;
    target_position = 2;
    sleeve = 1;
    add_reflector = 1;
    geo = get_Proto_geometry(1,1,skimmer,target_position,sleeve,add_reflector);    
    plot(Zeval,Reval,'ko')
    plot(f.z,f.r,'b','linewidth',2)
    plot(f.z(end),f.r(end),'rx','markersize',12)
end

%--------------------------------------------------------------------------
% Example 4: 
% Check if a field line intersects a component inside the vessel before 
% striking the targe
%--------------------------------------------------------------------------
Reval = 0.07; Zeval = 2.2;
dz_want = 0.01;
L = Ztarg - Zeval;
nsteps = round(abs(L/dz_want));
dz = L/nsteps;
f = follow_fieldlines_rzphi_dz(bfield,Reval,Zeval,0,dz,nsteps);
fprintf('\nExample 4\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,f.r(end),Ztarg)

skimmer = 1;
target_position = 2;
sleeve = 1;
add_reflector = 1;
geo = get_Proto_geometry(1,1,skimmer,target_position,sleeve,add_reflector);
plot(Zeval,Reval,'ko')
plot(f.z,f.r,'b','linewidth',2)
plot(f.z(end),f.r(end),'rx','markersize',12)

allin = all(inpolygon(f.z,f.r,geo.vessel_clip_z,geo.vessel_clip_r));
if allin
    fprintf('No intersections detected\n')
else
    fprintf('Found intersection(s)\n')
    ins = inpolygon(f.z,f.r,geo.vessel_clip_z,geo.vessel_clip_r);
    int1 = find(ins==0,1,'first') - 1;
    fprintf('Last point before intersection (R,Z) = (%6.3f,%6.3f) [m]\n',f.r(int1),f.z(int1))
    plot(f.z(int1),f.r(int1),'mx','markersize',12,'linewidth',3)
end


%--------------------------------------------------------------------------
% Example 5: 
% Make a 2D contour plot and highlight the contour=field line from a point
%--------------------------------------------------------------------------
nr = 100; nz = 500;
Z1d = linspace(0.0,5,nz);
R1d = linspace(0,0.1,nr);
% R1d = linspace(0,0.175,nr);
[R2d,Z2d] = meshgrid(R1d,Z1d);
psi2d = calc_psi_mpex(coil,current,R2d,Z2d);

% Reval = 0.045; Zeval = 1.75;
Reval = 0.0628  ; Zeval = 1.75;
psi_eval = calc_psi_mpex(coil,current,Reval,Zeval);

skimmer = 1;
target_position = 2;
sleeve = 1;
add_reflector = 1;
geo = get_Proto_geometry(1,1,skimmer,target_position,sleeve,add_reflector);
plot(Zeval,Reval,'ko')
contour(Z1d,R1d,psi2d.',logspace(-5.5,log10(max(max(psi2d))),20))
contour(Z1d,R1d,psi2d.',[1,1]*psi_eval,'k','linewidth',3)

