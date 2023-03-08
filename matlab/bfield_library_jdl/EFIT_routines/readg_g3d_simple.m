function g = readg_g3d_simple(filename)
% Simple reading of gfile
%
% Note: 
%   F = R*Btor
%  Jtor(Amp/m2) = R*P'(psi) + FF'(psi)/R
% JDL

%% Parse inputs
if ~isfile(filename)
    error('File not found: %s',filename)
end

fid = fopen(filename,'r');
if fid == -1
    error('Could not open file: %s\n',filename)
end

%% Line 1
% This should follow the formatting, but there are many files that do not.
% All we really need here are mw and mh, so there are a few attempts to get
% this from files that do not follow the original format.
% format (6a8,3i4)

line = fgetl(fid);
if length(line) < 60
    % strict formatting not followed
    % Try getting two integers from the end
    twoIntsRev = sscanf(fliplr(strtrim(line)),'%d %d',2);
    g.mh = str2double(fliplr(num2str(twoIntsRev(1))));
    g.mw = str2double(fliplr(num2str(twoIntsRev(2))));
else
    % Try formatted read
    g.ecase = sscanf(line(1:8),'%8s');
    sscanf(line(8*6+1:8*6+3),'%4i');   % time
    g.mw = sscanf(line(53:56),'%4i');  % Number of horizontal R grid points
    g.mh = sscanf(line(57:60),'%4i');  % Number of vertical Z grid points
end
% one last try if above failed
if (length(g.mw) ~= 1) || (length(g.mh) ~= 1) % some gfiles seem to have this shifted by a space or two
    temp = sscanf(line(50:end),'%d %d %d');
    g.mw = temp(2);
    g.mh = temp(3);
end
g.line1 = line;

%% Line 2
% Format 5e16.9
line = fgetl(fid);
data = sscanf(line,'%f%f%f%f%f',5);
g.xdim = data(1);   % rdim. Horizontal dim in meter of comp. box
g.zdim = data(2);   % zdim  Vertical dim in meter of comp. box
g.rzero = data(3);  % rcentr R in meter of toroidal magnetic field bcentr
g.rgrid1 = data(4); % rleft Min R in meter of comp box
g.zmid = data(5);   % zmid  Z in center of grid in m

%% Line 3
% Format 5e16.9
line = fgetl(fid);
data = sscanf(line,'%f%f%f%f%f',5);
g.rmaxis = data(1); % R of magnetic axis in meter
g.zmaxis = data(2); % Z of magnetic axis in meter
g.ssimag = data(3); % simag Poloidal flux at mag. axis in Weber/rad
g.ssibry = data(4); % sibry Poloidal flux at the plasma boundary in Weber/rad
g.bcentr = data(5); % Vacuum toroidal mag field in Tesla at RCENTR

%% Line 4
% Some files have empty lines for 4 and 5
% Format 5e16.9
test = fgetl(fid);
if ~isempty(test)
    line = sscanf(test,'%f',5);
    g.cpasma = line(1); % current Plasma current in Ampere
    % remaining data is redundant or meaningless (simag,xdum,rmaxis,xdum)
else
    g.cpasma = [];
end

%% Line 5
% Format 5e16.9
test = fgetl(fid);
if ~isempty(test)
    sscanf(test,'%f',5);
    % all data is redundant or meaningless (zmaxis,xdum,sibry,xdum,xdum)
end

%%
g.fpol = fscanf(fid,'%f',g.mw).';  % Poloidal current function in m-T. F = R*Bt
g.pres = fscanf(fid,'%f',g.mw).';  % Plasma pressure in Nt/m^2 on uniform flux grid
g.ffprim = fscanf(fid,'%f',g.mw).'; % FF'(psi) in (mT)^2/(Wb/rad) on uniform flux grid
g.pprime = fscanf(fid,'%f',g.mw).'; % P'(psi) in (Nt/m^2)/(Wb/rad) on uniform flux grid
g.psirz = fscanf(fid,'%f',[g.mw,g.mh]);  % Poloidal flux in Weber/rad on the rectangular grid points
g.qpsi = fscanf(fid,'%f',g.mw).'; % q values on uniform flux grid from axis to boundary

if any(g.qpsi == 0)
    g.qpsi = [];
end
if ~isempty(g.qpsi)  % file ends here in some cases
    %
    line = fscanf(fid,'%i',2);
    % Format 2i5
    g.nbdry = line(1);  % nbbbs number of boundary pts
    g.limitr = line(2); % number of limiter pts
    %
    g.bdry = fscanf(fid,'%f',[2,g.nbdry]); % rbbs, zbbs
    %
    g.lim = fscanf(fid,'%f',[2,g.limitr]); % rlim, zlim
    if max(g.lim(:)) > 100
        fprintf('Warning: lim seems to be in [cm], converting to [m]\n')
        g.lim = g.lim./100;
    end
else
    g.nbdry = 0;
    g.limitr = 0;
    g.bdry = [];
    g.lim = [];
end
fclose(fid);
