function g = readg_g3d(filename,force_rewrite,quiet)
% Simple reading of gfile with postprocessing for psi and bfield interpolation

if nargin < 2
    force_rewrite = 0;
end
if nargin < 3
    quiet = 0;
end

% Version identifies changes to routine that may require re-creating .mat file
version_ = 9;
% v3: Changed outputs to row vectors and added bicub inverse
% v4: Added aminor and nGW
% v5: More robust processing of incomplete eqdsk files, remove spline fits
% v6: removing old bicub matrix
% v7: Add tri calculation
% v8: Add better shape metrics
% v9: Add Bpol_est estimate from Ip filament at magnetic axis

if isempty(filename)
    g = [];
    return;
end

fname_mat = [filename,'.mat'];
iwarn = 0;
[action,ierr] = check_file_exist_and_new(filename,fname_mat,iwarn,0);
if ierr
    error('Could not find gfile or cached .mat file: %s',filename)
end

if force_rewrite
    action = 'raw';
end

if strcmp(action,'mat')
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
else
    if ~quiet
        disp([' >>>> Reading gfile ',filename])
    end
end

try
    g = readg_g3d_simple(filename);
catch ME
    rawFileAccessError = strncmp(ME.message,'File not found:',15) || ...
        strncmp(ME.message,'Could not open file:',20);
    if rawFileAccessError && ~force_rewrite && exist(fname_mat,'file') == 2
        if ~quiet
            fprintf('Info: Unable to read raw gfile, using .mat cache: %s\n',fname_mat)
            fprintf('      Raw read error: %s\n',ME.message)
        end
        S = load(fname_mat);
        g = S.g;
        g.filename = filename;
        return
    elseif rawFileAccessError && ~force_rewrite
        if ~quiet
            fprintf('Info: Unable to read raw gfile, returning []: %s\n',filename)
            fprintf('      Raw read error: %s\n',ME.message)
        end
        g = [];
        return
    end
    rethrow(ME)
end

%% Add shaping parameters and other grid details
g = postprocess_gfile(g);


%% Try parsing filename for shot and time
g.filename = filename;
try
    [~,f2,f3] = fileparts(filename);
    g.shot = sscanf(f2(2:end),'%d');
    g.gfilename = [f2,f3];
catch
end

g.version = version_;

if ~quiet
    disp('>>> Saving .mat version of gfile')
end
try
    save(fname_mat,'g');
catch
    if ~quiet
        fprintf('Info: Was unable to save .mat version of gfile: %s\n',fname_mat)
    end
end

end
