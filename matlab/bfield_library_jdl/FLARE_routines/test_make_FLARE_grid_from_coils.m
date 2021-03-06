clearvars;

RECALC_MAT = 1;



if RECALC_MAT
    if 0
        
        
        out_dir = 'C:\Work\FLARE\coils';
        % out_file = 'Bgrid_vac_22kA_mimic_45x41x32_v2'; % 22kA mimic  -- corrected IS1
        % out_file = 'Bgrid_vac_22kA_mimic_90x82x65_v2'; % 22kA mimic  -- corrected IS1
        out_file = 'Bgrid_vac_narrow_mirror_90x82x65'; % Narrow mirror
        out_file = fullfile(out_dir,out_file);
        
        coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
        coil = load_vmec_coils_file(coils_file);
        winding_array = [108,108,108,108,108,36,36,8,8];
        taper_norm = [1, 1.02, 1.08, 0.97, 0.88, 0.15, -0.15, 0, 0].*winding_array; Inorm = 12556;  % Narrow mirror
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, 0.0150, -0.0150]; Inorm = 1.341e6;%  0kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1950, -0.1050, 0.0150, -0.0150]; Inorm = 1.354e6; % 11kA mimic
%         taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1700, -0.1300, 0.0150, -0.0150]; Inorm = 1.367e6; % 22kA mimic  -- corrected IS1
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1450, -0.1550, 0.0150, -0.0150]; Inorm = 1.380e6; % 32kA mimic
        % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, 0.0150, -0.0150]; Inorm = 1.393e6; % 43kA mimic
        
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
        
    elseif 1

        xdr_path = 'C:\Work\Stellarator\W7X EMC3 modeling\BGRID\';

%         fname_xdr_dump = 'fieldn_mfbe181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_022.dat';
%         out_file = 'Bgrid_xdr_mfbe_22kA_90x82x65';
        
        fname_xdr_dump = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_022.dat';
        out_file = 'Bgrid_xdr_altern_22kA_90x82x65';

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

WRITE_DIR = 2;

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