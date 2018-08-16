function [num_lines,fdata] = num_lines_file(file_name)
    % Return number of lines in file, optionally file data as cell array 
    % with one line per element.
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
    if nargout > 1
       fdata = f{1}(:); 
    end

