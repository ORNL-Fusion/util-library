function wout=load_wout(wout_file,display_info)

if nargin < 2
    display_info = 0;
end

% % wout_file = 'C:\Work\Pellets\wout_vmec_pcnet.nc';

ncid = netcdf.open(wout_file,'NC_NOWRITE');
[ndim0,nvar0,natt0,nunlim0] = netcdf.inq(ncid);


for i = 0:ndim0-1
    [dimname,dimlength] = netcdf.inqDim(ncid,i);
    if display_info
        fprintf('Dimension %i: %20s, length: %i\n',i,dimname,dimlength)
    end
end
if display_info
    disp('-------------------------------------------')
end
% xtype: 1= byte, 2 = character, 3 = short integer, 4 = integer, 5 = real, 6 = double
for i = 0:nvar0-1
    [varname,xtype,dimid,natt] = netcdf.inqVar(ncid,i);
    varnames{i+1} = varname;
    if display_info
        fprintf('Variable %i: %20s, dim: ',i,varname)
        for j = 1:length(dimid)
            fprintf('%i     ',dimid(j))
        end
        fprintf('\n')
    end
end

varid = netcdf.inqVarID(ncid,'lasym__logical__');
wout.asym = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'extcur');
wout.extcur = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'version_');
wout.version = netcdf.getVar(ncid,varid);
wout.xm = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'xm');
wout.xm = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'rmax_surf');
wout.rmax_surf = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'rmin_surf');
wout.rmin_surf = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'zmax_surf');
wout.zmax_surf = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'xn');
wout.xn = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'mnmax');
wout.mnmax = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'mnmax_nyq');
wout.mnmax_nyq = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'xn_nyq');
wout.xn_nyq = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'xm_nyq');
wout.xm_nyq = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'ns');
wout.ns = double(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'presf');
wout.presf = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'iotaf');
wout.iotaf = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'iotas');
wout.iotas = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'phi');
wout.phi = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'rmnc');
wout.rmnc = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'zmns');
wout.zmns = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'lmns');
wout.lmns = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'gmnc');
wout.gmnc = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'bsupumnc');
wout.bsupumnc = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'bsupvmnc');
wout.bsupvmnc = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'bsubumnc');
wout.bsubumnc = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'bsubvmnc');
wout.bsubvmnc = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'bsubsmns');
wout.bsubsmns = netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'nfp');
wout.nfp = double(netcdf.getVar(ncid,varid));

if version > 8.5
    varid = netcdf.inqVarID(ncid,'equif');
    wout.equif = netcdf.getVar(ncid,varid);
end

if wout.asym
    varid = netcdf.inqVarID(ncid,'rmns');
    wout.rmns = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'zmnc');
    wout.zmnc = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'gmns');
    wout.gmns = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'lmnc');
    wout.lmnc = netcdf.getVar(ncid,varid);    
    varid = netcdf.inqVarID(ncid,'bsupumns');
    wout.bsupumns = netcdf.getVar(ncid,varid);    
    varid = netcdf.inqVarID(ncid,'bsupvmns');
    wout.bsupvmns = netcdf.getVar(ncid,varid); 
    varid = netcdf.inqVarID(ncid,'bsubumns');
    wout.bsubumns = netcdf.getVar(ncid,varid); 
    varid = netcdf.inqVarID(ncid,'bsubvmns');
    wout.bsubvmns = netcdf.getVar(ncid,varid); 
    varid = netcdf.inqVarID(ncid,'bsubsmnc');
    wout.bsubsmnc = netcdf.getVar(ncid,varid);     
end

wout.s = wout.phi./wout.phi(end);
wout.shalf = wout.s - 0.5*(wout.s(2)-wout.s(1));
% dphids = phi(end);

wout.splrmn = spline(wout.s,wout.rmnc);
wout.splzmn = spline(wout.s,wout.zmns);
wout.spllmn = spline(wout.s,wout.lmns);
wout.splgmn = spline(wout.s,wout.gmnc);
wout.splbsupu = spline(wout.shalf,wout.bsupumnc);
wout.splbsupv = spline(wout.shalf,wout.bsupvmnc);
wout.splbsubu = spline(wout.shalf,wout.bsubumnc);
wout.splbsubv = spline(wout.shalf,wout.bsubvmnc);
wout.splbsubs = spline(wout.shalf,wout.bsubsmns);
wout.spliotas = spline(wout.shalf,wout.iotas);
    
if wout.asym
    wout.splrmns = spline(wout.s,wout.rmns);
    wout.splzmnc = spline(wout.s,wout.zmnc);
    wout.spllmnc = spline(wout.s,wout.lmnc);
    wout.splgmns = spline(wout.s,wout.gmns);
    wout.splbsupus = spline(wout.shalf,wout.bsupumns);
    wout.splbsupvs = spline(wout.shalf,wout.bsupvmns);
    wout.splbsubus = spline(wout.shalf,wout.bsubumns);
    wout.splbsubvs = spline(wout.shalf,wout.bsubvmns);
    wout.splbsubsc = spline(wout.shalf,wout.bsubsmnc);
end


netcdf.close(ncid);
