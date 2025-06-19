clearvars;

this_dir = pwd; 

% Typically there is a data file that is slow to read, which is then saved 
% as a file with the same name with ".mat" appended, which can be loaded
% more quickly on subsequent reads. However, if the original file has 
% been updated, then we want to reread it and replace the .mat file.
% 

fname_raw     = fullfile(this_dir,'mydata');
fname_mat = strcat(fname_raw,'.mat');

iwarn = 1;
purge = 0;

M = magic(5);

% Test 1: 
%   Neither file exists
fprintf('\nTest 1, neither file exists.\n')
if isfile(fname_raw), delete(fname_raw); end
if isfile(fname_mat), delete(fname_mat); end
[action,ierr] = check_file_exist_and_new(fname_raw,fname_mat,iwarn,purge);
fprintf('Test result || action: %s, ierr: %d\n\n',action, ierr);

% Test 2: 
% Only file "mydata" exists
fprintf('Test 2, only raw exists.\n')
if isfile(fname_raw), delete(fname_raw); end
if isfile(fname_mat), delete(fname_mat); end
dlmwrite(fname_raw,M);
[action,ierr] = check_file_exist_and_new(fname_raw,fname_mat,iwarn,purge);
fprintf('Test result || action: %s, ierr: %d\n\n',action, ierr);

% Test 3: 
% Only file "mydata.mat" exists
fprintf('Test 3, only mat exists.\n')
if isfile(fname_raw), delete(fname_raw); end
if isfile(fname_mat), delete(fname_mat); end
save(fname_mat,'M');
[action,ierr] = check_file_exist_and_new(fname_raw,fname_mat,iwarn,purge);
fprintf('Test result || action: %s, ierr: %d\n\n',action, ierr);
% Test 4: 
% Both exist and "mydata" is newer
fprintf('Test 4, both exist and raw newer.\n')
if isfile(fname_raw), delete(fname_raw); end
if isfile(fname_mat), delete(fname_mat); end
save(fname_mat,'M');
pause(1)
dlmwrite(fname_raw,M);
[action,ierr] = check_file_exist_and_new(fname_raw,fname_mat,iwarn,purge);
fprintf('Test result || action: %s, ierr: %d\n\n',action, ierr);
% Test 5: 
% Both exist and "mydata.mat" is newer
fprintf('Test 5, both exist and mat newer.\n')
if isfile(fname_raw), delete(fname_raw); end
if isfile(fname_mat), delete(fname_mat); end
dlmwrite(fname_raw,M);
pause(1)
save(fname_mat,'M');
[action,ierr] = check_file_exist_and_new(fname_raw,fname_mat,iwarn,purge);
fprintf('Test result || action: %s, ierr: %d\n\n',action, ierr);
if isfile(fname_raw), delete(fname_raw); end
if isfile(fname_mat), delete(fname_mat); end
