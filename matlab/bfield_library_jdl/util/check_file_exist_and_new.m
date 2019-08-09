function [read_raw,read_mat,ierr] = check_file_exist_and_new(fname_raw,fname_mat,iwarn,purge)
% 
% Typically there is a data file that is slow to read, which is then saved 
% as a file with the same name with ".mat" appended, which can be loaded
% more quickly on subsequent reads. However, if the original file has 
% been updated, then we want to reread it and replace the .mat file.
% 

if nargout < 3
    error('Routine has been updated, check returns')
end

read_mat = [];
read_raw = [];
ierr = 1;

if exist(fname_raw,'file') == 2
    raw_exists = 1;
else
    if iwarn
        fprintf('Did not find raw output file: %s\n',fname_raw)
    end    
    raw_exists = 0;
end

if exist(fname_mat,'file') == 2
    if purge == 1
        delete(fname_mat);
        mat_exists = 0;
    else
        mat_exists = 1;
    end
else
    mat_exists = 0;
end

if ~raw_exists && ~mat_exists
    ierr = 1;
    read_raw = 0;
    read_mat = 0;    
    if iwarn
        fprintf('Did not find raw or mat output file: %s, %s\n',fname_raw,fname_mat)
    end
    return;
elseif raw_exists && ~mat_exists
    ierr = 0;
    read_raw = 1;
    read_mat = 0;    
    return;
elseif ~raw_exists && mat_exists
    ierr = 0;
    read_raw = 0;
    read_mat = 1;    
    if iwarn
        fprintf('No raw file but found mat version: %s, %s\n',fname_raw,fname_mat)
    end    
    return;
end

% Both files exist, check timestamps
finfo_raw = dir(fname_raw);
finfo_mat = dir(fname_mat);
if finfo_mat.datenum > finfo_raw.datenum
    ierr = 0;
    read_raw = 0;
    read_mat = 1;
else
    ierr = 0;
    read_raw = 1;
    read_mat = 0;    
end
