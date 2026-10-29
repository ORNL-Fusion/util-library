function g = readg_g3d(filename,force_rewrite,quiet)
% Simple reading of gfile with postprocessing for psi and bfield interpolation
% g.psiUnitsDiagnostic reports an advisory Wb versus Wb/rad comparison.

if nargin < 2
    force_rewrite = 0;
end
if nargin < 3
    quiet = 0;
end

% Version identifies changes to routine that may require re-creating .mat file
version_ = 10;
% v3: Changed outputs to row vectors and added bicub inverse
% v4: Added aminor and nGW
% v5: More robust processing of incomplete eqdsk files, remove spline fits
% v6: removing old bicub matrix
% v7: Add tri calculation
% v8: Add better shape metrics
% v9: Add Bpol_est estimate from Ip filament at magnetic axis
% v10: Add advisory gfile psi-unit diagnostic

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
            g = refresh_fpol_pchip(g);
            warn_if_gfile_psi_may_be_wb(g);
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
        g = refresh_fpol_pchip(g);
        warn_if_gfile_psi_may_be_wb(g);
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
warn_if_gfile_psi_may_be_wb(g);

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

function warn_if_gfile_psi_may_be_wb(g)
% Report a substantial 2*pi discrepancy on every read, including cached reads.
if ~isfield(g,'psiUnitsDiagnostic') || ~g.psiUnitsDiagnostic.available
    return
end
Diagnostic = g.psiUnitsDiagnostic;
if ~strcmp(Diagnostic.inferredUnits,'Wb')
    return
end
scoreWbPerRad = abs(log(Diagnostic.bpolWbPerRad/Diagnostic.bpolEstimate));
scoreWb = abs(log(Diagnostic.bpolWb/Diagnostic.bpolEstimate));
if scoreWbPerRad-scoreWb < log(2)
    return
end
fprintf(2,'WARNING: gfile %s: psirz may be in Wb rather than Wb/rad.\n',g.filename);
fprintf(2,'  |Ip|-based Bpol estimate at the outboard LCFS = %.6g T.\n',Diagnostic.bpolEstimate);
fprintf(2,'  Gfile gradient gives Bpol = %.6g T as Wb/rad, or %.6g T as Wb.\n',Diagnostic.bpolWbPerRad,Diagnostic.bpolWb);
fprintf(2,'  No flux conversion was applied; derived B fields assume Wb/rad.\n');
fprintf(2,'  If Wb is confirmed, call g = correct_gfile_wb_to_wb_per_radian(g).\n');
end

function g = refresh_fpol_pchip(g)
% Upgrade an in-memory v10 cache without rebuilding the flux grid.
g.fpol_pp = pchip(g.pn,g.fpol);
if isfield(g,'fpol_coeffs')
    g = rmfield(g,'fpol_coeffs');
end
fprintf(2,'  If Wb is confirmed, call g = correct_gfile_wb_to_wb_per_radian(g);\n');
end
