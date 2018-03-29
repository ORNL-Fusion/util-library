function num_lines = num_lines_file(file_name)
    if isstring(file_name) || ischar(file_name)
        if exist(file_name,'file') ~= 2
            error('File does not exist: %s\n',file_name)
        end
         fid = fopen(file_name,'r');
    else
        fid = file_name;
    end
    f = textscan(fid,'%s','delimiter','\n','whitespace','');
    num_lines = length(f{1});
    if isstring(file_name)
        fclose(fid);
    else
        frewind(fid);
    end

