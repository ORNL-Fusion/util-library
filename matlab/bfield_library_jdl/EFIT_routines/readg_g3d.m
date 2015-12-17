function g = readg_g3d(filename)

fname_mat = [filename,'.mat'];
if exist(fname_mat,'file') ~= 2
    disp([' >>>> Reading gfile ',filename])
else
    disp([' >>>> Reading .mat version of gfile ',fname_mat])    
    load(fname_mat);
    return;    
end

fid = fopen(filename,'r');

line = fgetl(fid); g.ecase = sscanf(line(1:8),'%8s'); sscanf(line(8*6+1:8*6+3),'%4i'); g.mw = sscanf(line(53:56),'%4i'); g.mh = sscanf(line(57:60),'%4i');
line = fscanf(fid,'%f%f%f%f%f',5); g.xdim = line(1); g.zdim = line(2); g.rzero = line(3); g.rgrid1 = line(4); g.zmid = line(5);
line = fscanf(fid,'%f%f%f%f%f',5); g.rmaxis = line(1); g.zmaxis = line(2); g.ssimag = line(3); g.ssibry = line(4); g.bcentr = line(5);
line = fscanf(fid,'%f',5); g.cpasma = line(1);
fscanf(fid,'%f',5);
g.fpol = fscanf(fid,'%f',g.mw);
g.pres = fscanf(fid,'%f',g.mw);
g.ffprim = fscanf(fid,'%f',g.mw);
g.pprime = fscanf(fid,'%f',g.mw);
g.psirz = fscanf(fid,'%f',[g.mw,g.mh]);
g.qpsi = fscanf(fid,'%f',g.mw);
line = fscanf(fid,'%i',2); g.nbdry = line(1); g.limitr = line(2);
g.bdry = fscanf(fid,'%f',[2,g.nbdry]);
g.lim = fscanf(fid,'%f',[2,g.limitr]);
fclose(fid);

%
% Postprocessing
% 
g.dR = g.xdim/(g.mw-1);
g.dZ = g.zdim/(g.mh-1);

for i = 0:g.mw-1
    g.r(i+1) = g.rgrid1 + g.dR*i;
    g.pn(i+1,1) = i/(g.mw-1);
end

for i = 0:g.mh-1
    g.z(i+1) = g.zmid - 0.5*g.zdim + g.dZ*i;
end

g.ip_sign = -sign(g.cpasma);

g.bicub_coeffs = get_psi_bicub_coeffs(g);
g.fpol_coeffs = polyfit(g.pn,g.fpol,7);

disp(' >> Saving .mat version of gfile')
save(fname_mat,'g');
