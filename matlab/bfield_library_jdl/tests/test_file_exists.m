clearvars;

this_dir = pwd; 

% Typically there is a data file that is slow to read, which is then saved 
% as a file with the same name with ".mat" appended, which can be loaded
% more quickly on subsequent reads. However, if the original file has 
% been updated, then we want to reread it and replace the .mat file.
% 

fname     = fullfile(this_dir,'mydata');
fname_mat = strcat(fname,'.mat');

iwarn = 1;
purge = 0;

M = magic(5);

% Test 1: 
%   Neither file exists
fprintf('\nTest 1, neither file exists.\n')
if isfile(fname), delete(fname); end
if isfile(fname_mat), delete(fname_mat); end
[read_raw,read_mat,ierr] = check_file_exist_and_new(fname,fname_mat,iwarn,purge);
fprintf('Test 1 result || read_raw: %d, read_mat: %d, ierr: %d\n\n',read_raw,read_mat, ierr);

% Test 2: 
% Only file "mydata" exists
fprintf('Test 2, only raw exists.\n')
if isfile(fname), delete(fname); end
if isfile(fname_mat), delete(fname_mat); end
dlmwrite(fname,M);
[read_raw,read_mat,ierr] = check_file_exist_and_new(fname,fname_mat,iwarn,purge);
fprintf('Test 2 result || read_raw: %d, read_mat: %d, ierr: %d\n\n',read_raw,read_mat, ierr);

% Test 3: 
% Only file "mydata.mat" exists
fprintf('Test 3, only mat exists.\n')
if isfile(fname), delete(fname); end
if isfile(fname_mat), delete(fname_mat); end
save(fname_mat,'M');
[read_raw,read_mat,ierr] = check_file_exist_and_new(fname,fname_mat,iwarn,purge);
fprintf('Test 3 result || read_raw: %d, read_mat: %d, ierr: %d\n\n',read_raw,read_mat, ierr);

% Test 4: 
% Both exist and "mydata" is newer
fprintf('Test 4, both exist and raw newer.\n')
if isfile(fname), delete(fname); end
if isfile(fname_mat), delete(fname_mat); end
save(fname_mat,'M');
pause(1)
dlmwrite(fname,M);
[read_raw,read_mat,ierr] = check_file_exist_and_new(fname,fname_mat,iwarn,purge);
fprintf('Test 4 result || read_raw: %d, read_mat: %d, ierr: %d\n\n',read_raw,read_mat, ierr);

% Test 5: 
% Both exist and "mydata.mat" is newer
fprintf('Test 5, both exist and mat newer.\n')
if isfile(fname), delete(fname); end
if isfile(fname_mat), delete(fname_mat); end
dlmwrite(fname,M);
pause(1)
save(fname_mat,'M');
[read_raw,read_mat,ierr] = check_file_exist_and_new(fname,fname_mat,iwarn,purge);
fprintf('Test 5 result || read_raw: %d, read_mat: %d, ierr: %d\n\n',read_raw,read_mat, ierr);

if isfile(fname), delete(fname); end
if isfile(fname_mat), delete(fname_mat); end
