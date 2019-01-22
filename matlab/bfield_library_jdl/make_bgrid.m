clearvars;

RECALC_MAT = 1;



if RECALC_MAT
    if 1
        out_dir = 'C:\Work\Stellarator\W7X EMC3 modeling\BGRID';

        

        coils_file = 'C:\Work_archive\RUN_ARCHIVE\VMEC_RUNS\W7X_mark_coils_and_input\coils.w7x';
        coil = load_vmec_coils_file(coils_file);
        winding_array = [108,108,108,108,108,36,36,8,8];
%         taper = [12022	11897	12148	13399	13524	8219	-3005  2500 -2500]; out_file = 'Bgrid_vac_0kA_mimic_90x82x65';  % 0kA mimic, run in OP1.2a EES+252
%         taper = [12129	12002	12255	13519	13645	7391	-3980  2500 -2500]; out_file = 'Bgrid_vac_11kA_mimic_90x82x65'; % 11kA mimic, run in OP1.2a EFS+252
%         taper = [12243	12116	12371	13646	13774	6504	-4973  2500 -2500]; out_file = 'Bgrid_vac_22kA_mimic_90x82x65'; % 22kA mimic, run in OP1.2a EGS+252
        taper = [12243	12116	12371	13646	13774	7754       -3723  2500 -2500]; out_file = 'Bgrid_vac_22kA_2pct_iota_red_mimic_90x82x65'; % 22kA mimic, run in OP1.2a EGS+252
%         taper = [12359	12230	12487	13775	13904	5600	-5987  2500 -2500]; out_file = 'Bgrid_vac_32kA_mimic_90x82x65'; % 32kA mimic, run in OP1.2a EGS001+252
%         taper = [12477	12347	12607	13907	14037	4679	-7019  2500 -2500]; out_file = 'Bgrid_vac_43kA_mimic_90x82x65'; % 43kA mimic, run in OP1.2a FHS+252
%           taper = [12477	12347	12607	13907	14037	5879	-5819  2500 -2500]; out_file = 'Bgrid_vac_43kA_2pct_iota_red_mimic_90x82x65'; % 43kA mimic, run in OP1.2a FHS+252
%          taper = [12477	12347	12607	13907	14037	4679	-7019  1250 -1250]; out_file = 'Bgrid_vac_43kA_+1250_mimic_90x82x65'; % 43kA mimic, run in OP1.2a FHS+252
%         taper = [12477	12347	12607	13907	14037	4679	-7019  0 0]; out_file = 'Bgrid_vac_43kA_0_mimic_90x82x65'; % 43kA mimic, run in OP1.2a FHS+252
%         taper = [12477	12347	12607	13907	14037	4679	-7019  -1250 1250]; out_file = 'Bgrid_vac_43kA_-1250_mimic_90x82x65'; % 43kA mimic, run in OP1.2a FHS+252
%         taper = [12477	12347	12607	13907	14037	4679	-7019  -2500 2500]; out_file = 'Bgrid_vac_43kA_-2500_mimic_90x82x65'; % 43kA mimic, run in OP1.2a FHS+252

%         taper = [1, 1, 1, 1, 1, 0, 0, 0, 0]*1e6; out_file = 'Bgrid_w7x_standard_1MA_90x82x65.mat'; % 
%         taper = [1, 1, 1, 1, 1, -0.18, -0.18, 0, 0]*1e6; out_file = 'Bgrid_w7x_low_iota_-180_1MA_90x82x65'; %        
%         taper = [1, 1, 1, 1, 1, 0.030, 0.030, 0, 0]*1e6; out_file = 'Bgrid_w7x_low_iota_30_1MA_90x82x65'; %
%         taper = [1, 1, 1, 1, 1, 0.105, 0.105, 0, 0]*1e6; out_file = 'Bgrid_w7x_low_iota_105_1MA_90x82x65'; %
%         taper = [1, 1, 1, 1, 1, 0.360, 0.360, 0, 0]*1e6; out_file = 'Bgrid_w7x_low_iota_360_1MA_90x82x65'; %
%         taper = [1, 1, 1, 1, 1, 0.480, 0.480, 0, 0]*1e6; out_file = 'Bgrid_w7x_low_iota_480_1MA_90x82x65'; %
%         taper = [1, 1, 1, 1, 1, 0.750, 0.750, 0, 0]*1e6; out_file = 'Bgrid_w7x_low_iota_750_1MA_90x82x65'; %
        coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
        out_file = fullfile(out_dir,out_file);
        
        bfield.type = 'just_coils';
        bfield.coil = coil.coil;
        bfield.current = coil.current;
        bfield.nsym = coil.num_periods;
        
    elseif 0
        
        
        out_dir = 'C:\Work\FLARE\coils';        
        % out_file = 'Bgrid_vac_22kA_mimic_45x41x32_v2'; % 22kA mimic  -- corrected IS1
        % out_file = 'Bgrid_vac_22kA_mimic_90x82x65_v2'; % 22kA mimic  -- corrected IS1
%         out_file = 'Bgrid_vac_narrow_mirror_90x82x65'; % Narrow mirror
        out_file = 'Bgrid_vac_43kA_mimic_+2500_90x82x65';
        out_file = fullfile(out_dir,out_file);
        
        coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
        coil = load_vmec_coils_file(coils_file);
        winding_array = [108,108,108,108,108,36,36,8,8];
%         taper_norm = [1, 1.02, 1.08, 0.97, 0.88, 0.15, -0.15, 0, 0].*winding_array; Inorm = 12556;  % Narrow mirror
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, 0.0150, -0.0150]; Inorm = 1.341e6;%  0kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1950, -0.1050, 0.0150, -0.0150]; Inorm = 1.354e6; % 11kA mimic
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1700, -0.1300, 0.0150, -0.0150]; Inorm = 1.367e6; % 22kA mimic  -- corrected IS1
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1450, -0.1550, 0.0150, -0.0150]; Inorm = 1.380e6; % 32kA mimic
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, 0.0150, -0.0150]; Inorm = 1.393e6; % 43kA mimic
        
        
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, 0.01436, -0.01436]; Inorm = 1.393e6; % 43kA mimic +2500
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, 0.007178750, -0.007178750]; Inorm = 1.393e6; % 43kA mimic +1250
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, +0.007178750, 0.007178750]; Inorm = 1.393e6; % 43kA mimic -1250
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, 0, 0]; Inorm = 1.393e6; % 43kA mimic +0
        
        
        %     taper_norm = [1,1,1,1,1,1,1,1]; Inorm = 1; out_file = 'Bgrid_coil1_1Amp.txt';
        taper = Inorm*taper_norm./winding_array;
        
        coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
        
        bfield.type = 'just_coils';
        bfield.coil = coil.coil;
        bfield.current = coil.current;
        bfield.nsym = coil.num_periods;
    elseif 0
        out_dir = 'C:\Work\Stellarator\W7X EMC3 modeling\VMEC_RUNS\W7X_2016\OP2_SE_ref\22kA\';
        out_file = 'Bgrid_bmw_90x82x65';
%         out_file = 'Bgrid_fromA_bmw_90x82x65';
        out_file = fullfile(out_dir,out_file);
        
        bmw_fname = fullfile(out_dir,'bmw_result.nc');
        bmw = load_bmw_file(bmw_fname);
        
        
        if 0
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
        
        
        bfield.nsym = bmw.nsym;
        
    elseif 0

        xdr_path = 'C:\Work\Stellarator\W7X EMC3 modeling\BGRID\';

%         fname_xdr_dump = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_000.dat';        
%         out_file = 'Bgrid_xdr_altern_OP2_00kA_90x82x65';
%         fname_xdr_dump = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_011.dat';        
%         out_file = 'Bgrid_xdr_altern_OP2_11kA_90x82x65';
%         fname_xdr_dump = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_022.dat';        
%         out_file = 'Bgrid_xdr_altern_OP2_22kA_90x82x65';
%         fname_xdr_dump = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_032.dat';        
%         out_file = 'Bgrid_xdr_altern_OP2_32kA_90x82x65';
        fname_xdr_dump = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_043.dat';        
        out_file = 'Bgrid_xdr_altern_OP2_43kA_90x82x65';        
        
%         fname_xdr_dump = 'fieldn_mfbe181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_022.dat';
%         out_file = 'Bgrid_xdr_mfbe_22kA_90x82x65';
        
%         fname_xdr_dump = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_022.dat';
%         out_file = 'Bgrid_xdr_altern_22kA_90x82x65';


        out_dir = xdr_path;

        
        out_file = fullfile(out_dir,out_file);        

        xdr = read_xdr_dump_file(xdr_path,fname_xdr_dump);
        
        bfield.type = 'xdr';
        bfield.nsym = xdr.nperio;
        bfield.xdr = xdr;
        
        
    else
        error('do something')
    end
    
%     Rmin = 4.25;
%     Rmax = 6.45;
%     Zmin = -1.2;
%     Zmax = 1.2;
    Rmin = 4.3;
    Rmax = 6.4;
    Zmin = -1.2;
    Zmax = 1.2;    
    
    nz = 45*2;  % Seems to work ok
    nr = 41*2; 
    nsym = 5;
    nphi = 33*2; % This is number of planes including 0 and 2*pi/5.  One fewer will be written to the data files (not the .mat file though)
    
    if exist(out_file,'file') == 2
        inp = input('The output file already exists, are you sure you want to overwrite [Y/N] ??');
        if ~strcmpi(inp,'y')
            error('Stopping')
        end
    end
    
    dr = (Rmax-Rmin)/(nr-1);
    dz = (Zmax-Zmin)/(nz-1);
    dphi = 2*pi/nsym/(nphi-1);
    fprintf('dR = %f\n',dr)
    fprintf('dZ = %f\n',dz)
    fprintf('dphi (deg) = %f\n',dphi*180/pi)
        
    r = linspace(Rmin,Rmax,nr);
    z = linspace(Zmin,Zmax,nz);
    p = linspace(0,2*pi/nsym,nphi);
    
    nowarn = 0;
    Br = zeros(nr,nz,nphi);
    Bz = zeros(nr,nz,nphi);
    Bphi = zeros(nr,nz,nphi);
    t0 = tic;
    for kt = 1:nphi
        fprintf('Working on kt %d of %d\n',kt,nphi)
        for jz = 1:nz
%             fprintf('    jz %d of %d\n',jz,nz)
            for ir = 1:nr
                [Bout,ierr] = bfield_general_rzphi(r(ir),z(jz),p(kt),bfield,nowarn);
%                 if ierr ~= 0
%                     a=1
%                 end
                Br(ir,jz,kt) = Bout.br;
                Bz(ir,jz,kt) = Bout.bz;
                Bphi(ir,jz,kt) = Bout.bphi;
            end
% aa=1;
        end
%         if toc(t0) > 10
    end
    
    
    Bgrid.Br = Br;
    Bgrid.Bz = Bz;
    Bgrid.Bphi = Bphi;
    Bgrid.nr = nr;
    Bgrid.nz = nz;
    Bgrid.nphi = nphi;
    Bgrid.nsym = nsym;
    Bgrid.Rmin = Rmin;
    Bgrid.Zmin = Zmin;
    Bgrid.Rmax = Rmax;
    Bgrid.Zmax = Zmax;
    Bgrid.R = r;
    Bgrid.Z = z;
    Bgrid.phi = p;
    Bgrid.dR = dr;
    Bgrid.dZ = dz;
    Bgrid.dphi = dphi;
    
    save([out_file,'.mat'],'Bgrid')
    
else
    load([out_file,'.mat'])
end

WRITE_DIR = 2; %Use 2!

fid = fopen([out_file,'_layout.dat'],'w');
fprintf(fid,'%d %d %d %d %f %f %f %f\n',Bgrid.nr,Bgrid.nz,Bgrid.nphi-1,Bgrid.nsym,Bgrid.Rmin,Bgrid.Rmax,Bgrid.Zmin,Bgrid.Zmax);
fclose(fid);

fid = fopen([out_file,'_r.dat'],'w');
fprintf('Writing Br\n')
if WRITE_DIR == 1
    write_loop1(fid,Bgrid.Br,Bgrid.nr,Bgrid.nz,Bgrid.nphi);
elseif WRITE_DIR == 2
    write_loop2(fid,Bgrid.Br,Bgrid.nr,Bgrid.nz,Bgrid.nphi);
end
fclose(fid);

fid = fopen([out_file,'_z.dat'],'w');
fprintf('Writing Bz\n')
if WRITE_DIR == 1
    write_loop1(fid,Bgrid.Bz,Bgrid.nr,Bgrid.nz,Bgrid.nphi);
elseif WRITE_DIR == 2
    write_loop2(fid,Bgrid.Bz,Bgrid.nr,Bgrid.nz,Bgrid.nphi);
end
fclose(fid);

fid = fopen([out_file,'_phi.dat'],'w');
fprintf('Writing Bphi\n')
if WRITE_DIR == 1
    write_loop1(fid,Bgrid.Bphi,Bgrid.nr,Bgrid.nz,Bgrid.nphi);
elseif WRITE_DIR == 2
    write_loop2(fid,Bgrid.Bphi,Bgrid.nr,Bgrid.nz,Bgrid.nphi);
end
fclose(fid);

function write_loop1(fid,B3D,nr,nz,nphi)
for k = 1:nphi-1
    for j = 1:nz
        for i = 1:nr
            fprintf(fid,'%18.12e\n',B3D(i,j,k));
        end
    end
end
end


function write_loop2(fid,B3D,nr,nz,nphi)
for i = 1:nr
    for j = 1:nz
        for k = 1:nphi-1
            fprintf(fid,'%18.12e\n',B3D(i,j,k));
        end
    end
end
end