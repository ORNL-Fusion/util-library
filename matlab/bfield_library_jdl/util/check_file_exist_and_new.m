function [action,ierr] = check_file_exist_and_new(fname_raw,fname_mat,iwarn,purge)
% 
% Typically there is a data file that is slow to read, which is then saved 
% as a file with the same name with ".mat" appended, which can be loaded
% more quickly on subsequent reads. However, if the original file has 
% been updated, then we want to reread it and replace the .mat file.
% 
%  action   'raw'  : read slow-to-parse original file and refresh cache
%           'mat'  : load cached .mat file
%           'none' : nothing usable (caller must handle)
%  ierr   0 = OK, 1 = error (no usable file)
%
% iwarn : print diagnostic messages
% purge: if true delete the .mat file

narginchk(2,4);
if nargin < 3
    iwarn = 1;   
end
if nargin < 4
    purge = false; 
end

action = 'none';
ierr = 0;

raw_exists = isfile(fname_raw);

if purge && isfile(fname_mat)
    delete(fname_mat);
end
mat_exists = isfile(fname_mat);



if ~raw_exists && ~mat_exists
    ierr = 1;                         % neither file present
    if iwarn
        fprintf('Did not find raw or mat output file: %s  |  %s\n',fname_raw,fname_mat);
    end
    return
elseif raw_exists && ~mat_exists
    action = 'raw';                     % parse raw
    return
elseif ~raw_exists && mat_exists
    action = 'mat';                     % raw missing
    if iwarn
        fprintf('No raw file, using mat version: %s\n',fname_mat);
    end
    return
end

% Both files exist, check timestamps
raw_time = dir(fname_raw).datenum;
mat_time = dir(fname_mat).datenum;

if mat_time > raw_time
    action = 'mat';                     % cache newer than raw
else
    action = 'raw';                     % raw newer, refresh cache
end