function mgrid = load_mgrid_file(fname,display_netcdf_info,display_mgrid_info,no_binfo)

if nargin < 2
    display_netcdf_info = 0;
end
if nargin < 3
    display_mgrid_info = 0;
end
if nargin < 4
    no_binfo = 0;
end


ncid = netcdf.open(fname,'NC_NOWRITE');
[ndim0,nvar0,natt0,nunlim0] = netcdf.inq(ncid);


for i = 0:ndim0-1
    [dimname,dimlength] = netcdf.inqDim(ncid,i);
    if display_netcdf_info
        fprintf('Dimension %i: %20s, length: %i\n',i,dimname,dimlength)
    end
end
if display_netcdf_info
    disp('-------------------------------------------')
end
% xtype: 1= byte, 2 = character, 3 = short integer, 4 = integer, 5 = real, 6 = double
for i = 0:nvar0-1
    [varname,xtype,dimid,natt] = netcdf.inqVar(ncid,i);
    varnames{i+1} = varname;
    if display_netcdf_info
        fprintf('Variable %i: %20s, dim: ',i,varname)
        for j = 1:length(dimid)
            fprintf('%i     ',dimid(j))
        end
        fprintf('\n')
    end
end

% DIMENSIONS
dimid = netcdf.inqDimID(ncid,'external_coils'); [~,mgrid.external_coils]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'rad'); [~,mgrid.rad]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'zee'); [~,mgrid.zee]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'phi'); [~,mgrid.phi]= netcdf.inqDim(ncid,dimid);

% VARIABLES
varid = netcdf.inqVarID(ncid,'rmin'); mgrid.rmin = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'rmax'); mgrid.rmax = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'zmin'); mgrid.zmin = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'zmax'); mgrid.zmax = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'ir'); mgrid.ir = double(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'jz'); mgrid.jz = double(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'kp'); mgrid.kp = double(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'nfp'); mgrid.nfp = double(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'nextcur'); mgrid.nextcur = double(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'mgrid_mode'); mgrid.mgrid_mode = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'coil_group'); mgrid.coil_group = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'raw_coil_cur'); mgrid.raw_coil_cur = netcdf.getVar(ncid,varid);


% actual bfield data
if ~no_binfo
    for i = 1:mgrid.external_coils
        if i < 10
            cstring = strcat('00',num2str(i));
        elseif i < 100
            cstring = strcat('0',num2str(i));
        else
            error('Too many coils')
        end
        name = strcat('br_',cstring); varid = netcdf.inqVarID(ncid,name); mgrid.br{i} = netcdf.getVar(ncid,varid);
        name = strcat('bp_',cstring); varid = netcdf.inqVarID(ncid,name); mgrid.bp{i} = netcdf.getVar(ncid,varid);
        name = strcat('bz_',cstring); varid = netcdf.inqVarID(ncid,name); mgrid.bz{i} = netcdf.getVar(ncid,varid);
    end
end

% Grid
mgrid.r = linspace(mgrid.rmin,mgrid.rmax,mgrid.ir);
mgrid.z = linspace(mgrid.zmin,mgrid.zmax,mgrid.jz);
mgrid.phi_rad = linspace(0,2*pi/mgrid.nfp,mgrid.kp);


netcdf.close(ncid);

if display_mgrid_info
    fprintf('\n--------------------------------------------------\n')
    fprintf('Rmin,Rmax = [%f,%f]\n',mgrid.rmin,mgrid.rmax)
    fprintf('Zmin,Zmax = [%f,%f]\n',mgrid.zmin,mgrid.zmax)
    fprintf('ir,jz,kp = [%d,%d,%d]\n',mgrid.ir,mgrid.jz,mgrid.kp)
    fprintf('nfp,nextcur = [%d,%d]\n',mgrid.nfp,mgrid.nextcur)
    fprintf('mgrid mode: %s\n',mgrid.mgrid_mode)
    fprintf('--------------------------------------------------\n')
end