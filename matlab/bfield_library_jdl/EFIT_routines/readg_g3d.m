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

if ~isempty(g.qpsi)  % Do not postprocess truncated file
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

