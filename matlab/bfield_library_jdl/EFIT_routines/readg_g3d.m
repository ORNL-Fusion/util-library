function g = readg_g3d(filename,force_rewrite,quiet)
% Simple reading of gfile with postprocessing for psi and bfield interpolation

if nargin < 2
    force_rewrite = 0;
end
if nargin < 3
    quiet = 0;
end

% Version identifies changes to routine that may require re-creating .mat file
version_ = 8;
% v3: Changed outputs to row vectors and added bicub inverse
% v4: Added aminor and nGW
% v5: More robust processing of incomplete eqdsk files, remove spline fits
% v6: removing old bicub matrix
% v7: Add tri calculation
% v8: Add better shape metrics

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

%% Add shaping parameters and other grid details
g = postprocess_gfile(g);


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

end
