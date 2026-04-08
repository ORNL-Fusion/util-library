clearvars;

%%
%-----------------------------------------------------
% 0. Set up magnetic field -- MPEX parameters
%-----------------------------------------------------
config_name = 'D1-1';
verbose = 1;
% [Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson(config_name,verbose);
[Coil,windingCurrent,CoilGeometry,currentPerWinding] = build_MPEX_coils_jackson_hybrid(config_name,[],3,3,verbose);
[rcoil,zcoil] = get_coil_cross_sections(CoilGeometry);
bfield.coil = Coil; bfield.current = windingCurrent; bfield.type = 'MPEX';
Geo = get_MPEX_geometry;

%% Example 0: plot geometry
figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14);
plot_coil_cross_section(rcoil,zcoil,0,currentPerWinding);
colorbar;
plot(Geo.Vessel.z,Geo.Vessel.r,'r-','LineWidth',2)


%% Example 1: Calculate psi at a point R,Z
%--------------------------------------------------------------------------
% Reval = 0.063; Zeval = 1.75;
Reval = 0.045; Zeval = 1.75;
psi_eval = calc_psi_mpex(Coil,windingCurrent,Reval,Zeval);
fprintf('\nExample 1\n')
fprintf('At point (R,Z) = (%6.3f,%6.3f) [m], psi = %6.3e [Wb]\n',Reval,Zeval,psi_eval)



%% Example 2: 
% Map a point R,Z to the a line at constant axial position (e.g., target)
% by interpolating axial flux.
%--------------------------------------------------------------------------
Zmap = Geo.Target.z(1);  % Set axial position for mapping
Reval = 0.05; Zeval = 0;
Rmap = map_pt_to_target_mpex(Coil,windingCurrent,Zmap,Reval,Zeval);
fprintf('\nExample 2\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,Rmap,Zmap)



%% Example 3: 
% Map a point R,Z to a line at constant axial position (e.g., target)
% by following a field line
%--------------------------------------------------------------------------

MAKE_PLOT = 1;
dz_want = 0.01;
Reval = 0.05; Zeval = 0;
L = Geo.Target.z(1) - Zeval;
nsteps = round(abs(L/dz_want));
dz = L/nsteps;
f = follow_fieldlines_rzphi_dz(bfield,Reval,Zeval,0,dz,nsteps);
fprintf('\nExample 3\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,f.r(end),Zmap)

if MAKE_PLOT
    figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14);
    plot_coil_cross_section(rcoil,zcoil,0,currentPerWinding);
    plot(Geo.Vessel.z,Geo.Vessel.r,'r-','LineWidth',2)
    colorbar;
    plot(Zeval,Reval,'ko')
    plot(f.z,f.r,'b','linewidth',2)
    plot(f.z(end),f.r(end),'rx','markersize',12)
end




%% Example 4: 
% Check if a field line intersects a component inside the vessel before 
% striking the target
%--------------------------------------------------------------------------
Reval = 0.063; Zeval = 0;
dz_want = 0.01;
Ztarg = Geo.Target.z(1);
L = Ztarg - Zeval;
nsteps = round(abs(L/dz_want));
dz = L/nsteps;

f = follow_fieldlines_rzphi_dz(bfield,Reval,Zeval,0,dz,nsteps);
fprintf('\nExample 4\n')
fprintf('Point (R,Z) = (%6.3f,%6.3f) [m] maps to radius %8.5f [m] at Z = %6.3f [m]\n',Reval,Zeval,f.r(end),Ztarg)

figure; set(gcf,'color','w'); box on; grid on; hold on; set(gca,'fontsize',14);
plot_coil_cross_section(rcoil,zcoil,0,currentPerWinding);
plot(Geo.Vessel.z,Geo.Vessel.r,'r-','LineWidth',2)
colorbar;
plot(Zeval,Reval,'ko')
plot(f.z,f.r,'b','linewidth',2)
plot(f.z(end),f.r(end),'rx','markersize',12)


allin = all(inpolygon(f.z,f.r,Geo.Vessel.z,Geo.Vessel.z));
if allin
    fprintf('No intersections detected\n')
else
    fprintf('Found intersection(s)\n')
    ins = inpolygon(f.z,f.r,Geo.Vessel.z,Geo.Vessel.r);
    int1 = find(ins==0,1,'first') - 1;
    fprintf('Last point before intersection (R,Z) = (%6.3f,%6.3f) [m]\n',f.r(int1),f.z(int1))
    plot(f.z(int1),f.r(int1),'mx','markersize',12,'linewidth',3)
end



%% Example 5: 
% Make a 2D contour plot and highlight the contour=field line from a point
%--------------------------------------------------------------------------
nr = 100; nz = 500;
Z1d = linspace(-8,8,nz);
R1d = linspace(0,1.8,nr);
% R1d = linspace(0,0.175,nr);
[R2d,Z2d] = meshgrid(R1d,Z1d);
psi2d = calc_psi_mpex(Coil,windingCurrent,R2d,Z2d);

% Reval = 10e-2; Zeval = -2.85;
Reval = 4.75*0.0254/2; Zeval = -2743.20e-3 ;
psi_eval = calc_psi_mpex(Coil,windingCurrent,Reval,Zeval);

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
    patch(zcoil(i,:),rcoil(i,:),currentPerWinding(i))
end
colorbar; cc=get(gca,'clim'); set(gca,'clim',[0,cc(2)]);
xlabel('Z relative to ECH center (m)')
ylabel('|B|_{ax} (T)')

figure; hold on; box on; grid on;
plot(Z,Bout.bz)
for i = 1:size(rcoil,1)
    patch(zcoil(i,:),rcoil(i,:),currentPerWinding(i))
end
colorbar; cc=get(gca,'clim'); set(gca,'clim',[0,cc(2)]);
xlabel('Z relative to ECH center (m)')
ylabel('B^z_{ax} (T)')
