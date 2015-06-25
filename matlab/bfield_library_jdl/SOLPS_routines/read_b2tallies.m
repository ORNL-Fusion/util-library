function read_b2tallies(fname,display_netcdf_info)
if nargin < 2
    display_netcdf_info = 0;
end

display_netcdf_info = 1;
fname='C:\Work\SOLPS\b2tallies.nc';



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
dimid = netcdf.inqDimID(ncid,'vreg');    [~,b2.vreg]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'xreg');    [~,b2.xreg]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'yreg');    [~,b2.yreg]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'ns');      [~,b2.ns]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'speclng'); [~,b2.speclng]= netcdf.inqDim(ncid,dimid);
dimid = netcdf.inqDimID(ncid,'time');    [~,b2.time]= netcdf.inqDim(ncid,dimid);

% VARIABLES
varid = netcdf.inqVarID(ncid,'times'); b2.times = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'nereg'); b2.nereg = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'fnayreg'); b2.fnayreg = netcdf.getVar(ncid,varid);

netcdf.close(ncid);



