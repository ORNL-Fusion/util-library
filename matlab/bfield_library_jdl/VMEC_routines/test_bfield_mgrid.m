clearvars;



display_netcdf_info = 0;
display_mgrid_info = 1;
no_binfo = 0;
no_ainfo = 0;
nowarn = 0;

winding_array = [108,108,108,108,108,36,36,8,8];
taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.2200, -0.0800, -0.0150, 0.0150]; Inorm = 1.341e6;  %  0kA mimic
% taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1950, -0.1050, -0.0150, 0.0150]; Inorm = 1.354e6;  % 11kA mimic
% taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1700, -0.1300, -0.0150, 0.0150]; Inorm = 1.367e6;  % 22kA mimic
% taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1450, -0.1550, -0.0150, 0.0150]; Inorm = 1.380e6;  % 32kA mimic
% taper_norm = [0.9600, 0.9500, 0.9700, 1.0700, 1.0800, 0.1200, -0.1800, -0.0150, 0.0150]; Inorm = 1.393e6;  % 43kA mimic
%     taper_norm = [1, 1, 1, 1, 1, 1, 1, 1, 1]; Inorm = 1;  %  TEST
taper = Inorm*taper_norm./winding_array;

Rtest = 6;
Ztest = 0.1;
phitest =375*pi/180;

if 1
    coils_file = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\coils.w7x';
    coil = load_vmec_coils_file(coils_file);
    coil = set_w7x_current(coil,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
    bfield.type = 'just_coils';
    bfield.coil = coil.coil;
    bfield.current = coil.current;
    bfield.nsym = coil.num_periods;
    B = bfield_general_rzphi(Rtest,Ztest,phitest,bfield);
    fprintf('%20s %8.4e %8.4e %8.4e\n','B from B-S:',B.br,B.bz,B.bphi)
end

if 0
    
    tic;
    fname = 'C:\Work\Stellarator\W7X EMC3 modeling\Mark coils and input\mgrid_w7x.nc';
    mgrid = load_mgrid_file(fname,display_netcdf_info,display_mgrid_info,no_binfo,no_ainfo);
    toc
    

    
    mgrid = set_w7x_current(mgrid,taper); % taper = [I1,I2,I3,I4,I5,IA,IB,IS1,IS2];
    
    
    
    
    
    bfield.type = 'Bgrid';
    bfield.Bgrid = mgrid;
    bfield.nsym = mgrid.nsym;
    B = bfield_general_rzphi(Rtest,Ztest,phitest,bfield);
    fprintf('%20s %8.4e %8.4e %8.4e\n','B from mgrid:',B.br,B.bz,B.bphi)
    
    % make_Bgrid_from_mgrid(mgrid)
    
    
    
    
    
    
    
end

%----------------------------
if 0
    [Brcoeff,Bzcoeff,Bphicoeff,Bspline_info] = prepare_Bgrid_splines(mgrid);
    [Br,Bz,Bphi]=bfield_Bspline(Rtest,Ztest,phitest,Brcoeff,Bzcoeff,Bphicoeff,Bspline_info);
    fprintf('%20s %8.4e %8.4e %8.4e\n','B from Bspline:',Br,Bz,Bphi);
    
    [Arcoeff,Azcoeff,Aphicoeff,spline_info] = prepare_Agrid_splines(mgrid);
    [Br,Bz,Bphi]=bfield_Aspline(Rtest,Ztest,phitest,Arcoeff,Azcoeff,Aphicoeff,spline_info);
    fprintf('%20s %8.4e %8.4e %8.4e\n','B from Aspline:',Br,Bz,Bphi);
end


fname = 'C:\Work\Stellarator\W7X EMC3 modeling\mimic_configs\VAC\0kA\bmw_result.nc';
bmw = load_bmw_file(fname);
[Arcoeff,Azcoeff,Aphicoeff,spline_info] = prepare_Agrid_splines(bmw);
[Br,Bz,Bphi]=bfield_Aspline(Rtest,Ztest,phitest,Arcoeff,Azcoeff,Aphicoeff,spline_info);
fprintf('%20s %8.4e %8.4e %8.4e\n','B from Abmw:',Br,Bz,Bphi);

bfield.type = 'Aspline';
bfield.Arcoeff = Arcoeff;
bfield.Azcoeff = Azcoeff;
bfield.Aphicoeff = Aphicoeff;
bfield.spline_info = spline_info;
bfield.nsym = bmw.nsym;
B = bfield_general_rzphi(Rtest,Ztest,phitest,bfield);
fprintf('%20s %8.4e %8.4e %8.4e\n','B from Abmw gen:',B.br,B.bz,B.bphi)