clearvars;


if 0
    out_dir = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\';
    % out_file = 'Bgrid_90x82x66.mat';
    
    coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
    coil = load_vmec_coils_file(coils_file);
    winding_array = [108,108,108,108,108,36,36,8,8];
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, -0.0150, 0.0150]; Inorm = 1.341e6; out_file = 'Bgrid_w7x_0kA_mimic_180x164x132.mat'; %  0kA mimic
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1950, -0.1050, -0.0150, 0.0150]; Inorm = 1.354e6; out_file = 'Bgrid_w7x_11kA_mimic_180x164x132.mat';  % 11kA mimic
    taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1700, -0.1300, -0.0150, 0.0150]; Inorm = 1.367e6; out_file = 'Bgrid_w7x_22kA_mimic_180x164x132.mat'; % 22kA mimic
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1450, -0.1550, -0.0150, 0.0150]; Inorm = 1.380e6; out_file = 'Bgrid_w7x_32kA_mimic_180x164x132.mat'; % 32kA mimic
    % taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, -0.0150, 0.0150]; Inorm = 1.393e6; out_file = 'Bgrid_w7x_43kA_mimic_180x164x132.mat'; % 43kA mimic
    
    taper = Inorm*taper_norm./winding_array;
    
    coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
    
    bfield.type = 'just_coils';
    bfield.coil = coil.coil;
    bfield.current = coil.current;
    bfield.nsym = coil.num_periods;
elseif 1
    out_dir = 'C:\Work\Stellarator\W7X EMC3 modeling\mimic_configs\VAC\0kA\';
    out_file = 'Bgrid_bmw_180x164_132.mat';
    bmw_fname = fullfile(out_dir,'bmw_result.nc');
    bmw = load_bmw_file(bmw_fname);
    [Arcoeff,Azcoeff,Aphicoeff,spline_info] = prepare_Agrid_splines(bmw);
    
    bfield.type = 'Aspline';
    bfield.Arcoeff = Arcoeff;
    bfield.Azcoeff = Azcoeff;
    bfield.Aphicoeff = Aphicoeff;
    bfield.spline_info = spline_info;
    bfield.nsym = bmw.nsym;
    
else
    error('do something')
end



% Rmin = 4.25;
Rmin = 4.3;
Rmax = 6.4;
Zmin = -1.15;
Zmax = 1.15;


nz = 45*2*2;
nr = 41*2*2;
nsym = 5;
nphi = 33*2*2; % This is number of planes including 0 and 2*pi/5.  One fewer will be written to the file  <---- not doing this right now!!!


out_file = fullfile(out_dir,out_file);
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
for kt = 1:nphi
    fprintf('Working on kt %d of %d\n',kt,nphi)
    for jz = 1:nz
        fprintf('    jz %d of %d\n',jz,nz)
        for ir = 1:nr
            [Bout,ierr] = bfield_general_rzphi(r(ir),z(jz),p(kt),bfield,nowarn);
            Br(ir,jz,kt) = Bout.br;
            Bz(ir,jz,kt) = Bout.bz;
            Bphi(ir,jz,kt) = Bout.bphi;
        end
    end
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

save(out_file,'Bgrid')
% fid = fopen('C:\Work\Btest.dat');
% fprintf(fid,'%d %d %d %d %f %f %f %f\n',nr,nz,nphi,nsym,Rmin,Rmax,Zmin,Zmax);


            