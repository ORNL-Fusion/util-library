function bmw = load_bmw_file(fname)

fprintf('Reading bmw file\n')
ncid = netcdf.open(fname,'NC_NOWRITE');

% DIMENSIONS
dimid = netcdf.inqDimID(ncid,'r');   [~,bmw.nr  ] = netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'z');   [~,bmw.nz  ] = netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'phi'); [~,bmw.nphi] = netcdf.inqDim(ncid,dimid);

% VARIABLES
varid = netcdf.inqVarID(ncid,'rmin'); bmw.rmin = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'rmax'); bmw.rmax = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'zmin'); bmw.zmin = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'zmax'); bmw.zmax = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'nfp');  bmw.nfp = double(netcdf.getVar(ncid,varid));

i = 1;
name = 'ar_grid'; varid = netcdf.inqVarID(ncid,name); bmw.Ar{i}   = netcdf.getVar(ncid,varid);
name = 'ap_grid'; varid = netcdf.inqVarID(ncid,name); bmw.Aphi{i} = netcdf.getVar(ncid,varid);
name = 'az_grid'; varid = netcdf.inqVarID(ncid,name); bmw.Az{i}   = netcdf.getVar(ncid,varid);

name = 'br_grid'; varid = netcdf.inqVarID(ncid,name); bmw.Br{i}   = netcdf.getVar(ncid,varid);
name = 'bp_grid'; varid = netcdf.inqVarID(ncid,name); bmw.Bphi{i} = netcdf.getVar(ncid,varid);
name = 'bz_grid'; varid = netcdf.inqVarID(ncid,name); bmw.Bz{i}   = netcdf.getVar(ncid,varid);


netcdf.close(ncid);

% Grid
bmw.R = linspace(bmw.rmin,bmw.rmax,bmw.nr);
bmw.Z = linspace(bmw.zmin,bmw.zmax,bmw.nz);
bmw.dphi = 2*pi/bmw.nfp/bmw.nphi;
for i = 1:bmw.nphi
    bmw.phi(i) = (i-1)*bmw.dphi;
end
bmw.dR = bmw.R(2)- bmw.R(1);
bmw.dZ = bmw.Z(2)- bmw.Z(1);
bmw.nsym = bmw.nfp;
bmw.scale_factor = 1;
