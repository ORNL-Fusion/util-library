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

%% Find the namelist section
line=fgetl(fid);
while ~feof(fid) && ~contains(line,namelist,'IgnoreCase',true)
    line=fgetl(fid);
end
if feof(fid)
    disp(['ERROR: Namelist: ' namelist ' not found!']);
end

%% Now read the namelist
while ~feof(fid)
    line=strtrim(fgetl(fid));
    if isempty(line) % Skim empty lines
        break;
    end
    if ~strcmp(line(1),'!') % Skim comment lines
        break;
    end
end
total_line = [];
while ~strcmp(strtrim(line),'/')
    % Collect the namelist into one line so assignments can be separated by
    % equals signs, independent of the original line breaks.
    % Check for comments
    exdex=strfind(line,'!');
    if ~isempty(exdex)
        line(exdex(1):end) = [];
    end
    total_line=[total_line strtrim(line) ' '];
    line=fgetl(fid);
end
fclose(fid);

%% Cleanup
total_line=regexprep(total_line,'([\n|\r])','','ignorecase'); %Get Rid of CR
total_line=regexprep(total_line,' T ',' 1 ','ignorecase'); %Get Rid of T
total_line=regexprep(total_line,' F ',' 0 ','ignorecase'); %Get Rid of F
total_line = strrep(total_line,'=',' = ');
total_line=regexprep(total_line,'\.true\.',' 1 ','ignorecase'); %Get Rid of .true.
total_line=regexprep(total_line,'\.false\.',' 0 ','ignorecase'); %Get Rid of .false.

%%
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
        
        % Handle right side: data runs from this equals sign to the next
        % variable name.  Character data may contain words, so use the last
        % word before the next equals sign as the next variable name.
        
        % Special case for underscore not caught by isletter -- allow varnames to have underscore
        isUnderscore = false(1,length(this_right));
        isUnderscore(strfind(this_right,'_')) = 1;
        
        iTest = isletter(this_right) | isUnderscore;
        
        find_words = find(diff(iTest) == 1) + 1;
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
        % Expand Fortran repeat syntax, e.g. 3*0.0 -> 0.0,0.0,0.0.
        while any(strfind(this_data,'*'))
            thisDataTmp = this_data; 
            thisDataTmp(strfind(thisDataTmp,' ')) = ',';
            comdexStar = strfind(thisDataTmp,',');
        
        
            stardex=strfind(thisDataTmp,'*');
            stardex = stardex(1);
            nRepeat = thisDataTmp(comdexStar(find(comdexStar < stardex,1,'last'))+1:stardex-1);  % should be an integer
            repeatVal = thisDataTmp(stardex+1:comdexStar(find(comdexStar > stardex,1,'first'))-1);  % should be an real or logical, I guess
            
            newStr = strcat(repeatVal,',');
            for iRpt = 2:str2num(nRepeat)
                newStr(end+1:end+length(repeatVal)+1) = strcat(repeatVal,',');
            end
            this_data = strcat(this_data(1:comdexStar(find(comdexStar < stardex,1,'last'))-1),newStr,this_data(comdexStar(find(comdexStar > stardex,1,'first'))+1:end));
            
        end
        
    end
    
    % Should now have this_var, and this_data
    if DEBUG
        fprintf('get_data_from_namelist parsed %s = %s\n',this_var,this_data)
    end
    if strcmpi(var_want,this_var)
        if DEBUG
            fprintf('get_data_from_namelist found requested variable %s\n',var_want)
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


%% Return if nothing is found
if exist('found') == 0
    mydata = [];
    return;
end
% Handle found data
% if found.icount == 1

mydata = NaN(var_want_size);
% A namelist variable can be defined in more than one assignment.  Usually
% there is one full-array assignment, but indexed partial assignments like
% conpar(1,2,1)=... arrive as multiple found instances and are stitched into
% mydata here.
for i = 1:found.icount
    
    data_raw = regexprep(found.data{i},',',' ','ignorecase'); % Remove commas
    % Check for character data
    if any(strfind(data_raw,''''))
        data_raw = regexprep(data_raw,'''',' ');
        data = split(strip(data_raw));
        isChar = 1;
    else
        data = sscanf(data_raw,'%e'); % get numeric data
        isChar = 0;
    end

    if DEBUG
        fprintf('get_data_from_namelist working on %s assignment %d of %d\n',var_want,i,found.icount)
        fprintf('  requested size [%s], parsed %d values, nind = %d\n',num2str(var_want_size),length(data),found.nind(i))
        if found.nopar(i)
            fprintf('  detected assignment without explicit indices\n')
        else
            fprintf('  detected assignment with indices (%s)\n',found.par_inside{i})
        end
    end
    
    if prod(var_want_size) == length(data)
        if DEBUG
            fprintf('get_data_from_namelist assigning full %s, size [%s]\n',var_want,num2str(var_want_size))
        end
        mydata = data;
    else
        % If we got here could be a partial definition
        if length(data) > prod(var_want_size)
            if isChar
                if prod(var_want_size) == 1
                    % Stitch a single character string back together --
                    % allows there to be spaces in a string
                    mydata = strip(data_raw);
                    return;
                else
                    error('Too much data??')            
                end
            else
                error('Too much data??')            
            end
        else
            % partial def
            % check for colons
            coldex = strfind(found.par_inside{i},':');
            if ~isempty(coldex)
                if DEBUG
                    fprintf('get_data_from_namelist detected colon indexing for %s(%s)\n',var_want,found.par_inside{i})
                end
                error('need to handle colons')
            end
            if found.nopar(i)
                % Fortran namelists allow unindexed partial array assignments.
                % Start at the first array element and continue in column-major
                % order until the supplied data are exhausted.
                if numel(data) > numel(mydata)
                    error('get_data_from_namelist: too much unindexed partial data for %s',var_want)
                end
                if DEBUG
                    fprintf('get_data_from_namelist assigning unindexed partial %s from start, %d values\n',var_want,numel(data))
                end
                mydata(1:numel(data)) = data;
            elseif length(var_want_size) == found.nind(i)
                % should be easy case
                ind_data_raw = regexprep(found.par_inside{i},',',' ','ignorecase'); % Remove commas
                ind_data = sscanf(ind_data_raw,'%d');
                if DEBUG
                    fprintf('get_data_from_namelist partial %s(%s), %d values, offsets [%s]\n', ...
                        var_want,found.par_inside{i},numel(data),num2str(ind_offsets))
                end
                if ndims(mydata) == 1
                    ind_start = sub2ind(var_want_size,ind_data + ind_offsets);
                elseif ndims(mydata) == 2 %#ok<ISMAT>
                    ind_start = sub2ind(var_want_size,ind_data(1)+ind_offsets(1),ind_data(2)+ind_offsets(2));
                elseif ndims(mydata) == 3
                    ind_start = sub2ind(var_want_size,ind_data(1)+ind_offsets(1),ind_data(2)+ind_offsets(2),ind_data(3)+ind_offsets(3));
                elseif ndims(mydata) == 4
                    ind_start = sub2ind(var_want_size,ind_data(1)+ind_offsets(1),ind_data(2)+ind_offsets(2),ind_data(3)+ind_offsets(3),ind_data(4)+ind_offsets(4));                    
                else
                    error('Extend this to higher dimensions')
                end
                if DEBUG
                    fprintf('get_data_from_namelist writing %s linear indices %d:%d\n', ...
                        var_want,ind_start,ind_start + numel(data)-1)
                end
                
                mydata(ind_start:ind_start + numel(data)-1) = data;
            else
                % implied loop? -- should still work with 1d indices like above, just have to generalize sub2ind
                if DEBUG
                    fprintf('get_data_from_namelist cannot assign %s: expected %d indices, found %d\n', ...
                        var_want,length(var_want_size),found.nind(i))
                    fprintf('  raw data: %s\n',strtrim(data_raw))
                end
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
