function g = readg_g3d_simple(filename)
% Simple reading of gfile

if ~isfile(filename)
    error('File not found: %s',filename)
end

fid = fopen(filename,'r');
if fid == -1
    error(['Could not open file: ',filename])
end

% Line 1, part 1
line1 = fgetl(fid);
g.ecase = sscanf(line1(1:8),'%8s');
sscanf(line1(8*6+1:8*6+3),'%4i');   % time
g.mw = sscanf(line1(53:56),'%4i');  % Number of horizontal R grid points
g.mh = sscanf(line1(57:60),'%4i');  % Number of vertical Z grid points
if (length(g.mw) ~= 1) || (length(g.mh) ~= 1) % some gfiles seem to have this shifted by a space or two
    temp = sscanf(line1(50:end),'%d %d %d');
    g.mw = temp(2);
    g.mh = temp(3);
end
g.line1 = line1;

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
if any(g.qpsi == 0)
    g.qpsi = [];
end
if ~isempty(g.qpsi)  % file ends here in some cases
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
