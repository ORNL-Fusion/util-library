clearvars;

%-----------------------------------------------------
% Set up magnetic field -- MPEX parameters
%-----------------------------------------------------
% 2) Specify a 21 element array corresponding to each coil winding current
fil = define_MPEX_coil_filaments; % Get areas
% J_helicon = [4.917e7,4.917e7,2.92e6,2.92e6,3.917e7,3.917e7];
% J_ECH = [1.75e8,1.75e8,4.14e7,3.89e7];
% J_ICH = [5.6e7,5.41e7,5.378e7,5.543e7,4.275e7];
% J_transport = [6.372e7,2.37e7,5.62e7];
% J_target = [4.64e7,3.2e7,6.48e7];

% JE options 70 GHz and 105 GHz
J_helicon = [4.92e7,4.92e7,2.92e6,2.92e6,3.92e7,3.92e7];
J_ECH = [8.5e7,8.5e7,5.6e7,4.3e7]; 
J_ICH = [5.83e7,5.41e7,5.35e7,5.98e7,4.28e7]; 
J_transport = [6.37e7,2.37e7,5.62e7];
J_target = [4.64e7,3.07e7,6.48e7];

J_array = [J_helicon,J_ECH,J_ICH,J_transport,J_target];
current_in = J_array.*fil.area./fil.nwind;


% Build coils
verbose = 1;
[coil,current] = build_MPEX_coils_jackson(current_in,verbose);

[rcoil,zcoil] = get_MPEX_coil_cross_sections;
figure; hold on; box on;
for i = 1:size(rcoil,1)
    plot(zcoil(i,:),rcoil(i,:),'r')
    patch(zcoil(i,:),rcoil(i,:),J_array(i)./1e7)
end
colorbar; cc=get(gca,'clim'); set(gca,'clim',[0,cc(2)]);
% aasdf
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
if 0
Ztarg = 4.33415;  % Set axial position for mapping
% Ztarg = 2.68;  % Set axial position for mapping
Rmap = map_pt_to_target_mpex(coil,current,Ztarg,Reval,Zeval);
fprintf('\nExample 2\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,Rmap,Ztarg)
end
%--------------------------------------------------------------------------
% Example 3: 
% Map a point R,Z to a line at constant axial position (e.g., target)
% by following a field line
%--------------------------------------------------------------------------
MAKE_PLOT = 1;


bfield.coil = coil; bfield.current = current; bfield.type = 'MPEX';
if 0
dz_want = 0.01;
L = Ztarg - Zeval;
nsteps = round(abs(L/dz_want));
dz = L/nsteps;
f = follow_fieldlines_rzphi_dz(bfield,Reval,Zeval,0,dz,nsteps);
fprintf('\nExample 3\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,f.r(end),Ztarg)
end
% if MAKE_PLOT
%     skimmer = 1;
%     target_position = 2;
%     sleeve = 1;
%     add_reflector = 1;
%     geo = get_Proto_geometry(1,1,skimmer,target_position,sleeve,add_reflector);    
%     plot(Zeval,Reval,'ko')
%     plot(f.z,f.r,'b','linewidth',2)
%     plot(f.z(end),f.r(end),'rx','markersize',12)
% end
if 0
figure; hold on; box on; grid on;

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
end
% skimmer = 1;
% target_position = 2;
% sleeve = 1;
% add_reflector = 1;
% geo = get_Proto_geometry(1,1,skimmer,target_position,sleeve,add_reflector);
% plot(Zeval,Reval,'ko')
% plot(f.z,f.r,'b','linewidth',2)
% plot(f.z(end),f.r(end),'rx','markersize',12)
% figure; hold on; box on; grid on;

% allin = all(inpolygon(f.z,f.r,geo.vessel_clip_z,geo.vessel_clip_r));
% if allin
%     fprintf('No intersections detected\n')
% else
%     fprintf('Found intersection(s)\n')
%     ins = inpolygon(f.z,f.r,geo.vessel_clip_z,geo.vessel_clip_r);
%     int1 = find(ins==0,1,'first') - 1;
%     fprintf('Last point before intersection (R,Z) = (%6.3f,%6.3f) [m]\n',f.r(int1),f.z(int1))
%     plot(f.z(int1),f.r(int1),'mx','markersize',12,'linewidth',3)
% end


%--------------------------------------------------------------------------
% Example 5: 
% Make a 2D contour plot and highlight the contour=field line from a point
%--------------------------------------------------------------------------
nr = 100; nz = 500;
Z1d = linspace(-8,8,nz);
R1d = linspace(0,1.8,nr);
% R1d = linspace(0,0.175,nr);
[R2d,Z2d] = meshgrid(R1d,Z1d);
psi2d = calc_psi_mpex(coil,current,R2d,Z2d);

% Reval = 10e-2; Zeval = -2.85;
Reval = 4.75*0.0254/2; Zeval = -2743.20e-3 ;
psi_eval = calc_psi_mpex(coil,current,Reval,Zeval);

% skimmer = 1;
% target_position = 2;
% sleeve = 1;
% add_reflector = 1;
% geo = get_Proto_geometry(1,1,skimmer,target_position,sleeve,add_reflector);
plot(Zeval,Reval,'ko')
contour(Z1d,R1d,psi2d.',logspace(-5.5,log10(max(max(psi2d))),20))
aaa= contour(Z1d,R1d,psi2d.',[1,1]*psi_eval,'k','linewidth',3);

nz = 5000;

Z =  linspace(-6,8,nz);
R = 0.001*ones(size(Z));
P = 0*ones(size(Z));
[Bout,ierr] = bfield_general_rzphi(R,Z,P,bfield,1);

Bnorm = sqrt(Bout.bphi.^2 + Bout.br.^2 + Bout.bz.^2);
figure; hold on; box on; grid on;
plot(Z,Bnorm)
for i = 1:size(rcoil,1)
    patch(zcoil(i,:),rcoil(i,:),J_array(i)./1e7)
end
colorbar; cc=get(gca,'clim'); set(gca,'clim',[0,cc(2)]);
xlabel('Z relative to ECH center (m)')
ylabel('|B|_{ax} (T)')

figure; hold on; box on; grid on;
plot(Z,Bout.bz)
for i = 1:size(rcoil,1)
    patch(zcoil(i,:),rcoil(i,:),J_array(i)./1e7)
end
colorbar; cc=get(gca,'clim'); set(gca,'clim',[0,cc(2)]);
xlabel('Z relative to ECH center (m)')
ylabel('B^z_{ax} (T)')