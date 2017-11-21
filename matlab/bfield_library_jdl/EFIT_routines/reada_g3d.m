function a = reada_g3d(filename)
% filename = 'C:\Work\DIII-D\148712\a148712.04140_501';
fname_mat = [filename,'.mat'];
if exist(fname_mat,'file') ~= 2
    disp([' >>>> Reading afile ',filename])
else
    disp(['>>>> Reading .mat version of afile ',fname_mat])    
    S = load(fname_mat);
    a = S.a;
    return;    
end

fid = fopen(filename,'r');
if fid == -1
    error(['Could not open file: ',filename])
end
line = fgetl(fid); % datestr
line = fgetl(fid); tmp = sscanf(line,'%d%d'); a.shot = tmp(1); %tmp(2) = nt???
line = fgetl(fid); a.time_ms = sscanf(line,'%e');
line = fgetl(fid); % time, shape etc
line = fgetl(fid); tmp = sscanf(line,'%e%e%e%e'); a.tsaisq = tmp(1); a.rcencm = tmp(2); a.bcentr = tmp(3); a.pasmat = tmp(4); 
line = fgetl(fid); tmp = sscanf(line,'%e%e%e%e'); a.cpasma = tmp(1); a.rout   = tmp(2); a.zout   = tmp(3); a.aout   = tmp(4); 
line = fgetl(fid); tmp = sscanf(line,'%e%e%e%e'); a.eout   = tmp(1); a.doutu  = tmp(2); a.doutl  = tmp(3); a.vout   = tmp(4); 
line = fgetl(fid); tmp = sscanf(line,'%e%e%e%e'); a.rcurrt = tmp(1); a.zcurrt = tmp(2); a.qsta   = tmp(3); a.betat  = tmp(4); 
line = fgetl(fid); tmp = sscanf(line,'%e%e%e%e'); a.betap  = tmp(1); a.ali    = tmp(2); a.oleft  = tmp(3); a.oright = tmp(4); 

fclose(fid);

a.filename = filename;

disp('>>> Saving .mat version of afile')
save(fname_mat,'a');

