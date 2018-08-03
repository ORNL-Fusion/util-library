clearvars;

FIND_AXIS = 1;
% LOAD_POINCARE = 1;
BFIELD = -1;




switch BFIELD
    case -1
        poincare_out = 'low_iota_90_vac_poincare_0deg_50surf_100pt.mat';
        
        grid_file = 'C:\Work\Stellarator\W7X EMC3 modeling\BGRID\Bgrid_w7x_low_iota_90_1MA_90x82x65.mat';
        
%         grid_file = 'C:\Work\Stellarator\W7X EMC3 modeling\BGRID\Bgrid_vac_0kA_mimic_90x82x65.mat';
%             grid_file = 'C:\Work\Stellarator\W7X EMC3 modeling\BGRID\Bgrid_vac_43kA_mimic_90x82x65.mat';
        bdata = load(grid_file);
        bfield.type = 'Bgrid';
        bfield.Bgrid = bdata.Bgrid;
        bfield.nsym  = bdata.Bgrid.nsym;
        
        
    case 0
        poincare_out = 'delete_me_22kA_ref_coils_poincare_0deg_101surf_100pt.mat';
        
        coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
        coil = load_vmec_coils_file(coils_file);
        winding_array = [108,108,108,108,108,36,36,8,8];
        
        % First make poincare with vacuum coils, OP2 22kA
%         taper_norm = [12130, 12000, 12255, 13541, 13635, 9000, -2900, 0, 0].*winding_array; Inorm = 1; Rax = 5.915735; Zax = 0.000000; % 22kA OP2
        taper_norm = [1, 1.02, 1.08, 0.97, 0.88, 0.15, -0.15, 0, 0].*winding_array; Inorm = 12556; Rax = 5.931; Zax = 0.000000; % Narrow mirror
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, 0.0150, -0.0150]; Inorm = 1340931;%  0kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1950, -0.1050, 0.0150, -0.0150]; Inorm = 1353630; % 11kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1700, -0.1300, 0.0150, -0.0150]; Inorm = 1366546; % 22kA mimic  -- corrected IS1
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1450, -0.1550, 0.0150, -0.0150]; Inorm = 1379685; % 32kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, 0.0150, -0.0150]; Inorm = 1393054; % 43kA mimic
        
        taper = Inorm*taper_norm./winding_array;
        coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
        
        bfield.type = 'just_coils';
        bfield.coil = coil.coil;
        bfield.current = coil.current;
        bfield.nsym = coil.num_periods;
    case 1
        poincare_out = '22kA_ref_bmw_poincare_0deg_101surf_100pt.mat';
        
        if 0
            grid_file = 'C:\Work\Stellarator\W7X EMC3 modeling\VMEC_RUNS\W7X_2016\OP2_SE_ref\22kA\Bgrid_bmw_90x82x65.mat';
            %         grid_file = 'C:\Work\FLARE\coils\Bgrid_vac_22kA_mimic_90x82x65_v2.mat';
            bdata = load(grid_file);
            bfield.type = 'Bgrid';
            bfield.Bgrid = bdata.Bgrid;
            bfield.nsym  = bdata.Bgrid.nsym;
        else
            bmw_fname = 'C:\Work\Stellarator\W7X EMC3 modeling\VMEC_RUNS\W7X_2016\OP2_SE_ref\22kA\bmw_result.nc';
            bmw = load_bmw_file(bmw_fname);
        end
        if 1
            bmw.nphi = bmw.nphi + 1;
            bmw.phi(bmw.nphi) = 2*pi/bmw.nsym;
            bmw.Br{1}(:,:,bmw.nphi)   = bmw.Br{1}(:,:,1);
            bmw.Bz{1}(:,:,bmw.nphi)   = bmw.Bz{1}(:,:,1);
            bmw.Bphi{1}(:,:,bmw.nphi) = bmw.Bphi{1}(:,:,1);
            bfield.type = 'Bgrid';
            bfield.Bgrid = bmw;
            bfield.nsym  = bmw.nsym;            
        elseif 0
            [Arcoeff,Azcoeff,Aphicoeff,spline_info] = prepare_Agrid_splines(bmw);
            bfield.type       = 'Aspline';
            bfield.Arcoeff    = Arcoeff;
            bfield.Azcoeff    = Azcoeff;
            bfield.Aphicoeff  = Aphicoeff;
            bfield.spline_info = spline_info;
        else
            [Brcoeff,Bzcoeff,Bphicoeff,spline_info] = prepare_Bgrid_splines(bmw);
            bfield.type       = 'Bspline';
            bfield.Brcoeff    = Brcoeff;
            bfield.Bzcoeff    = Bzcoeff;
            bfield.Bphicoeff  = Bphicoeff;            
            bfield.spline_info = spline_info;
        end
        Rax = 5.915735; Zax = 0.000000; % 22kA OP2
        
        case 2
            grid_file = 'C:\Work\Stellarator\W7X EMC3 modeling\BGRID\Bgrid_xdr_altern_22kA_90x82x65.mat';
            temp = load(grid_file);
            bfield.Bgrid = temp.Bgrid;
            bfield.type = 'Bgrid';
            bfield.nsym = temp.Bgrid.nsym;

    otherwise
        error('bad BFIELD')
end

phistart = 0*pi/180; Rax_guess = 5.92; Zax_guess = 0;
% phistart = 36*pi/180; Rax_guess = 5.16; Zax_guess = 0;

if FIND_AXIS
    [Rax,Zax] = find_axis(bfield,phistart,Rax_guess,Zax_guess);
    [Bax] = bfield_general_rzphi(Rax,Zax,phistart,bfield);
    B0 = sqrt(Bax.br.^2 + Bax.bz.^2 + Bax.bphi.^2);
    fprintf('B0 = %f\n',B0)
end

% afsasdf

nsurf = 50;

if 0
    Rstart = 5.66; Rend = 5.60;
% Rstart = 6.1; Rend = 6.3;
    Zstart = 0; Zend = 0;
    Rstart = linspace(Rstart,Rend,nsurf);
    Zstart = linspace(Zstart,Zend,nsurf);
elseif 1
    % Define by distance from axis to point
    Rend = 5.255;
    Zend = 1.14;
    L = 0.0;  % m from ax to start surfs from
    Rstart=Rax-L/sqrt(1+((Zend-Zax)/(Rend-Rax))^2);
    Zstart=L/sqrt(1+((Rend-Rax)/(Zend-Zax))^2)-Zax;
    Rstart = linspace(Rstart,Rend,nsurf);
    Zstart = linspace(Zstart,Zend,nsurf);
end
    
poincare_out_dir = 'C:\Work\Stellarator\W7X EMC3 modeling\Poincare_data\';

poincare_save_name = fullfile(poincare_out_dir,poincare_out);

npoints_want = 100;
dphi = 1*pi/180;
sort_it = 1;

plot_settings.plotit = 1;
plot_settings.newfig = 1;
plot_settings.connect = 0;

if 1
    poinc = make_poincare(bfield,Rstart,Zstart,phistart,phistart,npoints_want,dphi,sort_it,Rax,Zax,plot_settings,poincare_save_name);
end
ves = load_W7X_vessel(0,0,phistart);
plot(ves.cut.r,ves.cut.z,'k')


asdfadf
run_path = 'C:\Work\EMC3\EMC3_runs\runs\W7X\mimic_22kA\VAC\block1_expand7\baserun';
run_info.run_path = run_path;
lim = load_all_limiter_files(run_info);
plot_emc3_plates_at_phi(lim,phistart);


asdfadf
wout_file = 'C:\Work\Stellarator\W7X EMC3 modeling\mimic_configs\VAC\0kA\wout_w7x.0990.1010.1115.1125.+0688.-0250.-0211.vacuum.0kA_mimic.v01.nc';
wout = load_wout(wout_file);
surf = plot_vmec_surfaces(wout,phistart,0);
plot(surf.rsurf(1:10:end,:).',surf.zsurf(1:10:end,:).','k','linewidth',2)