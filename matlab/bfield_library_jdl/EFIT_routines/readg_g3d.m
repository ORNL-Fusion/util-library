function g = readg_g3d(filename,force_rewrite)
% Simple reading of gfile with postprocessing for psi and bfield interpolation

if nargin < 2
    force_rewrite = 0;
end

% Version identifies changes to routine that may require re-creating .mat file
version_ = 3;
% v3: Changed outputs to row vectors and added bicub inverse

if isempty(filename)
    g = [];
    return;
end

fname_mat = [filename,'.mat'];
if exist(fname_mat,'file') ~= 2 || force_rewrite
    disp([' >>>> Reading gfile ',filename])
else
    disp(['>>>> Reading .mat version of gfile ',fname_mat])
    S = load(fname_mat);
    g = S.g;
    
    if isfield(g,'version')
        if g.version == version_
            return;
        end
    else
        g.version = 0;
    end
    fprintf(' Found version %i (or no version), expected %i, remaking.\n',g.version,version_)
end

fid = fopen(filename,'r');
if fid == -1
    error(['Could not open file: ',filename])
end

% Line 1, part 1
line = fgetl(fid);
g.ecase = sscanf(line(1:8),'%8s');
sscanf(line(8*6+1:8*6+3),'%4i');   % time
g.mw = sscanf(line(53:56),'%4i');  % Number of horizontal R grid points
g.mh = sscanf(line(57:60),'%4i');  % Number of vertical Z grid points
if (length(g.mw) ~= 1) || (length(g.mh) ~= 1) % some gfiles seem to have this shifted by a space or two
    temp = sscanf(line(50:end),'%d %d %d');
    g.mw = temp(2);
    g.mh = temp(3);
end

% Line 2
line = fgetl(fid);
data = sscanf(line,'%f%f%f%f%f',5);
g.xdim = data(1);
g.zdim = data(2);
g.rzero = data(3);
g.rgrid1 = data(4);
g.zmid = data(5);
% Line 3
line = fgetl(fid);
data = sscanf(line,'%f%f%f%f%f',5);
g.rmaxis = data(1);
g.zmaxis = data(2);
g.ssimag = data(3);
g.ssibry = data(4);
g.bcentr = data(5);
% Line 4
% Some files have empty lines for 4 and 5
test = fgetl(fid);
if ~isempty(test)
    line = sscanf(test,'%f',5);
    g.cpasma = line(1);
end
% Line 5
test = fgetl(fid);
if ~isempty(test)
    sscanf(test,'%f',5);
end
%
g.fpol = fscanf(fid,'%f',g.mw).';  % Poloidal current function in m-T = RBt
%
g.pres = fscanf(fid,'%f',g.mw).';  % Plasma pressure in Nt/m2
%
g.ffprim = fscanf(fid,'%f',g.mw).'; % FF'(psi) in (mT)%2/(Wb/rad)
%
g.pprime = fscanf(fid,'%f',g.mw).'; % P'(psi) in (Nt/m2)/(Wb/rad)
%
g.psirz = fscanf(fid,'%f',[g.mw,g.mh]);  % Poloidal flux in Weber/rad
%
g.qpsi = fscanf(fid,'%f',g.mw).';
if ~isempty(g.qpsi)  % file ends here
    %
    line = fscanf(fid,'%i',2);
    g.nbdry = line(1);
    g.limitr = line(2);
    %
    g.bdry = fscanf(fid,'%f',[2,g.nbdry]);
    %
    g.lim = fscanf(fid,'%f',[2,g.limitr]);
    if max(g.lim(:)) > 100
        fprintf('Warning: lim seems to be in [cm], converting to [m]\n')
        g.lim = g.lim./100;
    end
end
fclose(fid);

%
% Postprocessing
%
g.dR = g.xdim/(g.mw-1);
g.dZ = g.zdim/(g.mh-1);

g.r = zeros(1,g.mw);
g.z = zeros(1,g.mh);
g.pn = zeros(1,g.mw);
for i = 0:g.mw-1
    g.r(1,i+1) = g.rgrid1 + g.dR*i;
    g.pn(1,i+1) = i/(g.mw-1);
end

for i = 0:g.mh-1
    g.z(1,i+1) = g.zmid - 0.5*g.zdim + g.dZ*i;
end

if ~isempty(g.qpsi)  % file ends here
    g.ip_sign = -sign(g.cpasma);
    
    g.bicub_coeffs = get_psi_bicub_coeffs(g);
    g.bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g);
    g.fpol_coeffs = polyfit(g.pn,g.fpol,7);
    g.xk = dbsnak(g.mw,g.pn,3);
    g.fpol_coeffs_spline = dbsint(g.mw,g.pn,g.fpol,3,g.xk);
    g.filename = filename;
    
    %try to parse filename for shot and time
    [~,f2,f3] = fileparts(filename);
    g.shot = sscanf(f2(2:end),'%d');
    g.gfilename = [f2,f3];
end
g.version = version_;

disp('>>> Saving .mat version of gfile')
save(fname_mat,'g');

