function num_lines = num_lines_file(file_name)
if exist(file_name,'file') ~= 2
    error('File does not exist: %s\n',file_name)
end
fid = fopen(file_name,'r');
f = textscan(fid,'%s','delimiter','\n','whitespace','');
fclose(fid);
num_lines = length(f{1});
