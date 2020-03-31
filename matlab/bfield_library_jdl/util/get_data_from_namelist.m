function mydata = get_data_from_namelist(filename,namelist,var_want,var_want_size,ind_offsets)
% Still needs some extensions, but first pass
% 2019 JDL
if nargin < 5
    ind_offsets = zeros(1,length(var_want_size));
end

DEBUG = 0;

if ~isfile(filename)
    error('Could not find %s\n',filename);
end
fid=fopen(filename,'r');
% Find the namelist section
line=fgetl(fid);
while ~feof(fid) && ~contains(line,namelist)
    line=fgetl(fid);
end
if feof(fid)
    disp(['ERROR: Namelist: ' namelist ' not found!']);
end
% Now read the namelist
% Skim comment lines
while ~feof(fid)
    line=strtrim(fgetl(fid));
    if ~strcmp(line(1),'!')
        break;
    end
end
total_line = [];
while ~strcmp(strtrim(line),'/')
    % Check for comments
    exdex=strfind(line,'!');
    if ~isempty(exdex)
        line(exdex(1):end) = [];
    end
    total_line=[total_line strtrim(line) ' '];
    line=fgetl(fid);
end
fclose(fid);
% Cleanup
total_line=regexprep(total_line,'([\n|\r])','','ignorecase'); %Get Rid of CR
total_line=regexprep(total_line,' T ',' 1 ','ignorecase'); %Get Rid of T
total_line=regexprep(total_line,' F ',' 0 ','ignorecase'); %Get Rid of F
total_line=regexprep(total_line,' \.true\. ',' 1 ','ignorecase'); %Get Rid of .true.
total_line=regexprep(total_line,' \.false\. ',' 0 ','ignorecase'); %Get Rid of .false.
total_line = strrep(total_line,'=',' = ');

var_want = strtrim(var_want);
eqdex=strfind(total_line,'=');
istart = 1;
icount = 0;
for i = 1:length(eqdex)
    nind = 0;
    this_left  = total_line(istart:eqdex(i)-1);
    if i == length(eqdex)
        this_data = total_line(eqdex(i)+1:end);
    else
        this_right = total_line(eqdex(i)+1:eqdex(i+1)-1);
        
        % Handle right --> get data from eq to next varname
        find_words = find(diff(isletter(this_right)) == 1) + 1;
        if isempty(find_words)
            error('Could not find any letter strings in right side: %s',this_right)
        end
        % If there is only one word, it should be the next variable name
        % If there is more than one, then they should be character data
        %   in this case the next variable name should be the last word
        this_data = this_right(1:find_words(end)-1);
        istart = eqdex(i) + find_words(end);
    end
    
    % Handle left --> Look for varname
    % need to check for parentheses
    parleftdex = strfind(this_left,'(');
    parrightdex = strfind(this_left,')');
    if isempty(parleftdex) && isempty(parrightdex)
        % no parentheses
        this_var = lower(strtrim(this_left));
        nopar = 1;
        par_inside = [];
    elseif length(parleftdex) == 1 && length(parrightdex) == 1
        nopar = 0;
        par_inside = this_left(parleftdex+1:parrightdex-1);
        if isempty(par_inside)
            error('empty parentheses?')
        end
        comdex = strfind(par_inside,',');
        if isempty(comdex)
            nind = 1;
        else
            nind = length(comdex) + 1;
        end
        this_var = lower(strtrim(this_left(1:parleftdex-1)));
    else
        error('Do not understand parentheses structure')
    end
    
    if strfind(this_data,'*')
        error('need to handle this')
    end
    
    % Should now have this_var, and this_data
    if DEBUG
        fprintf('This var: %s\n',this_var)
        fprintf('This data: %s\n',this_data)
    end
    if strcmpi(var_want,this_var)
        if DEBUG
            fprintf('Found an instance of %s\n',var_want)
        end
        icount = icount + 1;
        found.eqdex(icount) = eqdex(i);
        found.nopar(icount) = nopar;
        found.data{icount} = this_data;
        found.nind(icount) = nind;
        found.par_inside{icount} = par_inside;
        found.icount = icount;
        
    end
    
end


% Handle found data
% if found.icount == 1

mydata = NaN(var_want_size);
for i = 1:found.icount
    
    data_raw = regexprep(found.data{i},',',' ','ignorecase'); % Remove commas
    % Check for character data
    if any(strfind(data_raw,''''))
        data_raw = regexprep(data_raw,'''',' ');
        data = split(strip(data_raw));
    else
        data = sscanf(data_raw,'%e'); % get numeric data
    end
    
    if prod(var_want_size) == length(data)
        mydata = data;
    else
        % If we got here could be a partial definition
        if length(data) > prod(var_want_size)
            error('Too much data??')
        else
            % partial def
            % check for colons
            coldex = strfind(par_inside,':');
            if ~isempty(coldex)
                error('need to handle colons')
            end
            if length(var_want_size) == found.nind(i)
                % should be easy case
                ind_data_raw = regexprep(found.par_inside{i},',',' ','ignorecase'); % Remove commas
                ind_data = sscanf(ind_data_raw,'%d');
                if ndims(mydata) == 1
                    ind_start = sub2ind(var_want_size,ind_data + ind_offsets);
                elseif ndims(mydata) == 2 %#ok<ISMAT>
                    ind_start = sub2ind(var_want_size,ind_data(1)+ind_offsets(1),ind_data(2)+ind_offsets(2));
                elseif ndims(mydata) == 3
                    ind_start = sub2ind(var_want_size,ind_data(1)+ind_offsets(1),ind_data(2)+ind_offsets(2),ind_data(3)+ind_offsets(3));
                else
                    error('more dims')
                end
                
                mydata(ind_start:ind_start + numel(data)-1) = data;
            else
                % implied loop? -- should still work with 1d indices like above, just have to generalize sub2ind
                error('extend this')
            end
            
        end
        
    end
end
% else
%     error('need to handle multi line')
%     % make sure to check for overwrites
% end
end