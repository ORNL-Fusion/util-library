function g = readg_g3d(filename,force_rewrite,quiet)
% Simple reading of gfile with postprocessing for psi and bfield interpolation

if nargin < 2
    force_rewrite = 0;
end
if nargin < 3
    quiet = 0;
end

% Version identifies changes to routine that may require re-creating .mat file
version_ = 6;
% v3: Changed outputs to row vectors and added bicub inverse
% v4: Added aminor and nGW
% v5: More robust processing of incomplete eqdsk files, remove spline fits
% v6: removing old bicub matrix

if isempty(filename)
    g = [];
    return;
end

fname_mat = [filename,'.mat'];
if exist(fname_mat,'file') ~= 2 || force_rewrite
    if ~quiet
        disp([' >>>> Reading gfile ',filename])
    end
else
    if ~quiet
        disp(['>>>> Reading .mat version of gfile ',fname_mat])
    end
    S = load(fname_mat);
    g = S.g;

    g.filename = filename;

    if isfield(g,'version')
        if g.version == version_
            return;
        end
    else
        g.version = 0;
    end
    fprintf(' Found version %i (or no version), expected %i, remaking.\n',g.version,version_)
end

g = readg_g3d_simple(filename);

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

g.aminor = [];
if ~isempty(g.bdry)
    g.aminor = (max(g.bdry(1,:))-min(g.bdry(1,:)))/2;    
end

g.nGW = [];
if ~isempty(g.cpasma) && ~isempty(g.aminor) 
    g.nGW = g.cpasma/1e6/pi/g.aminor^2*1e20;
end

if isempty(g.cpasma)
    g.ip_sign = 1;
else
    g.ip_sign = -sign(g.cpasma);
end

%% Flux interpolation coefficients
g.bicub_coeffs_inv = get_psi_bicub_coeffs_inv(g);

%% Fit F
g.fpol_coeffs = polyfit(g.pn,g.fpol,7);

%% Try parsing filename for shot and time
g.filename = filename;
try
    [~,f2,f3] = fileparts(filename);
    g.shot = sscanf(f2(2:end),'%d');
    g.gfilename = [f2,f3];
end

g.version = version_;

if ~quiet
    disp('>>> Saving .mat version of gfile')
end
save(fname_mat,'g');

