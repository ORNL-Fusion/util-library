clearvars;




if 1
    coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
    nfp_coil_plot = 5;
    coil = load_vmec_coils_file(coils_file);
    if 0
        query_vmec_coils_file(coil);
        plot_vmec_coils(coil,nfp_coil_plot);
    end
    
    
    winding_array = [108,108,108,108,108,36,36,8,8];
    taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, -0.0150, 0.0150]; Inorm = 1.341e6;  %  0kA mimic
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1950, -0.1050, -0.0150, 0.0150]; Inorm = 1.354e6;  % 11kA mimic
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1700, -0.1300, -0.0150, 0.0150]; Inorm = 1.367e6;  % 22kA mimic
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1450, -0.1550, -0.0150, 0.0150]; Inorm = 1.380e6;  % 32kA mimic
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, -0.0150, 0.0150]; Inorm = 1.393e6;  % 43kA mimic    
%     taper_norm = [1, 1, 1, 1, 1, 1, 1, 1, 1]; Inorm = 1;  %  TEST
    
    taper = Inorm*taper_norm./winding_array;
    
    coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
    bfield.type = 'just_coils';
    bfield.coil = coil.coil;
    bfield.current = coil.current;
    bfield.nsym = coil.num_periods;
else
    temp = load('C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\Bgrid_w7x_0kA_mimic_180x164x132.mat');
    Bgrid = temp.Bgrid;
    
    bfield.type = 'Bgrid';
    bfield.Bgrid = Bgrid;
    bfield.nsym = Bgrid.nsym;
end

if 1
    bounds.type = 'box';
    bounds.Rlims = [4.25,6.5];
    bounds.Zlims = [-1.3,1.3];
else
    ves = load_W7X_vessel(0,0);
    bounds.type = 'ves';
    bounds.ves = ves;
end
bfield.bounds = bounds;

phistart = 0*pi/180; Rax_guess = 5.92; Zax_guess = 0;
% phistart = 36*pi/180; Rax_guess = 5.16; Zax_guess = 0;


%
% Find axis and B0
%

[Rax,Zax] = find_axis(bfield,phistart,Rax_guess,Zax_guess);
[Bax] = bfield_general_rzphi(Rax,Zax,phistart,bfield);
B0 = sqrt(Bax.br.^2 + Bax.bz.^2 + Bax.bphi.^2);
fprintf('B0 = %f\n',B0)

%
% Poincare
%
nsurf = 3;

if 0
    Rstart = linspace(Rax+.275,Rax+1,nsurf);
    Zstart = Zax;
elseif 0
    % Define by distance from axis to point
    Rend = 5.3;
    Zend = 0.968;
    L = 0.2;  % m from ax to start surfs from
    Rstart=Rax-L/sqrt(1+((Zend-Zax)/(Rend-Rax))^2);
    Zstart=L/sqrt(1+((Rend-Rax)/(Zend-Zax))^2)-Zax;
    Rstart = linspace(Rstart,Rend,nsurf);
    Zstart = linspace(Zstart,Zend,nsurf);
else
    %     Rstart = 5.6;
    %     Zstart = 0;
    Rstart = linspace(5.9,5.95,3);
    Zstart = zeros(size(Rstart));
end

poincare_save_name = 'C:\Work\Stellarator\W7X EMC3 modeling\mimic_configs\VAC\0kA\poincare_0deg_3surf_25pt.mat';

npoints_want = 25;
dphi = 0.5*pi/180;
sort_it = 1;

plot_settings.plotit = 1;
plot_settings.newfig = 1;
plot_settings.connect = 0;

if 1
    poinc = make_poincare(bfield,Rstart,Zstart,phistart,phistart,npoints_want,dphi,sort_it,Rax,Zax,plot_settings,poincare_save_name);
end
ves = load_W7X_vessel(0,0,phistart);
plot(ves.cut.r,ves.cut.z,'k')


wout_file = 'C:\Work\Stellarator\W7X EMC3 modeling\mimic_configs\VAC\0kA\wout_w7x.0990.1010.1115.1125.+0688.-0250.-0211.vacuum.0kA_mimic.v01.nc';
wout = load_wout(wout_file);
surf = plot_vmec_surfaces(wout,phistart,0);
plot(surf.rsurf(1:10:end,:).',surf.zsurf(1:10:end,:).','k','linewidth',2)

% iuse = 16;
% R_surf = poinc.Rpoinc(:,iuse);
% Z_surf = poinc.Zpoinc(:,iuse);
% phi_surf = poinc.phi;
% phi_tor = calculate_toroidal_flux(bfield,R_surf,Z_surf,phi_surf);