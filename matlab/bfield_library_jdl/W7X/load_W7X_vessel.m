function ves = load_W7X_vessel(plotit,newfig,thisphi_rad)
% returns one FP of the W7X vessel in [m,rad]
% if thisphi_rad [rad] is provided then a single cut is returned interpolated
% at this angle

% vessel file assumed to be in cm, degrees, and is assumed to describe the
% first full field period including end points [0,72] degrees

mypath = fileparts(mfilename('fullpath'));  % Assumes vessel file is in the same dir as this file
fname = fullfile(mypath,'vessel.part');

if nargin < 1
    plotit = 0;
    newfig = 0;
end
if nargin == 1
    newfig = 1;
end
if nargin < 3
    thisphi_rad = [];
end

% Load vessel
fid = fopen(fname);

ves.label =fgetl(fid);
stuff = fscanf(fid,'%d %d %d %f %f %f %f',7);
ves.ntor = stuff(1);
ves.npol = stuff(2);
ves.nsym = stuff(3);
if ves.nsym ~= 5
    error('Double check file nsym ~= 5')
end
% rshift = stuff(4);
% zshift = stuff(5);

phipart = zeros(ves.ntor,ves.npol);
rpart   = zeros(ves.ntor,ves.npol);
zpart   = zeros(ves.ntor,ves.npol);
for i = 1:ves.ntor
    phitmp = fscanf(fid,'%f',1)*pi/180;
    for j = 1:ves.npol
        phipart(i,j) = phitmp;
        RZ = fscanf(fid,'%f %f %d',3);  % VESSEL FILE HAS 3 ENTRIES
        rpart(i,j) = RZ(1)/100;
        zpart(i,j) = RZ(2)/100;
    end
end
fclose(fid);

ves.r = rpart;
ves.z = zpart;
ves.x = ves.r.*cos(phipart);
ves.y = ves.r.*sin(phipart);
ves.phi = phipart;

if abs(ves.phi(1,1) - 0) > 1e-3 && abs(ves.phi(ves.ntor,1) - 2*pi/ves.nsym) > 1e-3
    error('Vessel should be defined from [0,2*pi/nsym]')
end

if thisphi_rad == 2*pi/ves.nsym    
    thisphi_fp1 = thisphi_rad;
else
    thisphi_fp1 = mod(thisphi_rad,2*pi/ves.nsym);
end

if ~isempty(thisphi_rad)
    ves.cut.x = zeros(1,ves.npol);
    ves.cut.y = zeros(1,ves.npol);
    ves.cut.z = zeros(1,ves.npol);
    ves.cut.r = zeros(1,ves.npol);
    ves.cut.phi = zeros(1,ves.npol);
    for j=1:ves.npol
        ves.cut.z(j) = interp1(ves.phi(:,j),ves.z(:,j),thisphi_fp1);
        ves.cut.r(j) = interp1(ves.phi(:,j),ves.r(:,j),thisphi_fp1);        
    end
    ves.cut.phi(:) = thisphi_rad;       
    ves.cut.x = ves.cut.r.*cos(thisphi_rad);
    ves.cut.y = ves.cut.r.*sin(thisphi_rad);   
end


if newfig == 1 
    figure;hold on;box on;
end
if plotit == 1
    plot3(ves.x,ves.y,ves.z,'b')
    plot3(ves.x.',ves.y.',ves.z.','b')
    if isfield(ves,'cut')
        plot3(ves.cut.x,ves.cut.y,ves.cut.z,'r.')
    end
end


