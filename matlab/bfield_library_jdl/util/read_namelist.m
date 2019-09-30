function data=read_namelist(filename,namelist)
%data=READ_NAMELIST(filename,namelist) Returns values from FORTRAN namelist
%   This function reads a FORTRAN input namelist and returns the values as
%   the fields of a structure.  Multidimensional arrays have their indicies
%   scaled and shifted to fit the matlab numbering scheme.
%
%   Note:  At this time the function does not support the colon modifier
%   and variables such as TEST(1:,5) = 1 2 3 4 will be overlooked.
%
%   Example:
%       data=read_namelist('input.test','INDATA');
%
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.1
%   Date:           1/26/11
%
% Extended by J.Lore
%   Arrays are assumed to start from index 1. If 0 is detected then arrays are shifted up by one

DEBUG = 1;

if DEBUG 
    fprintf('------------------------------------------\n')
    fprintf(' Working on initial read\n')
    fprintf('------------------------------------------\n')
end

% Open a text file
fid=fopen(filename,'r');
% Find the namelist section
line=fgetl(fid);
while ~feof(fid) && isempty(strfind(line,namelist)) 
    line=fgetl(fid);
end
if feof(fid)
    disp(['ERROR: Namelist: ' namelist ' not found!']);
    data=-1;
    return
end
data.temp=-1;
% Now read the namelist
while ~feof(fid)
    line=strtrim(fgetl(fid));
    if ~strcmp(line(1),'!')  % Skim comment lines
        break;
    end    
end
% total_line=line;
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
% Close the text file
fclose(fid);
% Now we extract every declaration
expr=['(([a-zA-Z0-9_]+)|([a-zA-Z0-9_]+\([\+\-0-9 ,]+\)))\s+='...
    '(((\s+(([0-9]+[*])?(\+|\-)?[0-9.]+E(\+|\-)?[0-9]+)+)+)|'...
    '(\s+\''.+\'')|'...
    '(\s+[T|F|\.true\.|\.false\.])|'...
    '((\s+([0-9]+[*])?(\+|\-)?[0-9.]+)+))'];
total_line=regexprep(total_line,'([\n|\r])','','ignorecase'); %Get Rid of CR
total_line=regexprep(total_line,' T ',' 1 ','ignorecase'); %Get Rid of T
total_line=regexprep(total_line,' F ',' 0 ','ignorecase'); %Get Rid of F
total_line=regexprep(total_line,' \.true\. ',' 1 ','ignorecase'); %Get Rid of .true.
total_line=regexprep(total_line,' \.false\. ',' 0 ','ignorecase'); %Get Rid of .false.
% Catch cases where poor formatting is used
total_line = strrep(total_line,'=',' = ');
% Also remove commas that are not inside of parentheses
eqdex=strfind(total_line,'=');
for i = 1:length(eqdex)
    % Look between the equals
    if i == length(eqdex)
        iend = length(total_line);
    else
        iend = eqdex(i+1)-1;
    end
    istart = eqdex(i)+1;
    part_check = total_line(istart:iend);
    % Check for commas
    comadex=strfind(part_check,',');
    if ~isempty(comadex) %#ok<*STREMP>
        % If there is a start parentheses, then shift string check
        par2dex=strfind(part_check,'(');
        if ~isempty(par2dex)
            iend = par2dex - 1;
            part_check = total_line(istart:iend);
            continue
        end
        % replace all commas
        total_line(istart:iend) = strrep(total_line(istart:iend),',',' ');
%         part_replace = strrep(part_check,',',' ');
%         total_line(istart+1:iend) = [];
%         total_line(istart) = '!';
%         total_line = strrep(total_line,'!',part_replace);       
    end
    
end

if DEBUG 
    fprintf('------------------------------------------\n')
    fprintf(' Working on regex and parsing\n')
    fprintf('------------------------------------------\n')
end
lines=regexpi(total_line,expr,'match');
% Now we parse each declaration
nlines=numel(lines);
for i=1:nlines
    line=handlestars(lines{i});
    if DEBUG
        fprintf('Working on parsed line %s\n',line);
    end
    
    
    eqdex=strfind(line,'=');
    par1dex=strfind(line,'(');
    par2dex=strfind(line,')');
    comadex=strfind(line,',');
    if isempty(comadex)
        fprintf('   >>> No commas in variable definition\n') % So 1D array assumed, with index taken to start from one
        if isempty(par1dex)
            fprintf('   >>> No parentheses in variable definition\n')
            name=lower(strtrim(sscanf(line(1:eqdex-1),'%s')));
            vals=sscanf(line(eqdex+1:numel(line)),'%g')';
            if isempty(vals)
                vals=strtrim(sscanf(line(eqdex+1:numel(line)),'%s'));
                vals=vals(2:numel(vals)-1); %Need to do this to handle ' in string
            end
            data.(name)=vals;
        else
            fprintf('   >>> Parentheses found in variable definition\n')  % So 1D array assumed with some indices
            name=lower(strtrim(sscanf(line(1:par1dex-1),'%s')));
            index=sscanf(line(par1dex+1:par2dex-1),'%g');
            vals=sscanf(line(eqdex+1:numel(line)),'%g')';
            if ~isfield(data,name)
                data.(name)=vals;
                data.([name '_nmlindex'])=index;
            else
                data.(name)=[data.(name) vals];
                data.([name '_nmlindex'])=[data.([name '_nmlindex']) index];
            end
        end
    else % Commas found     
        if DEBUG
            fprintf('   >>> Found commas in variable definition\n')
        end
        index_string='%g';
        for j=1:numel(comadex)
            index_string=[index_string ',%g'];
        end
        name=lower(strtrim(sscanf(line(1:par1dex-1),'%s')));
        fprintf('   >>> Variable name: %s\n',name);
        index=sscanf(line(par1dex+1:par2dex-1),index_string);
        vals=sscanf(line(eqdex+1:numel(line)),'%g')';
        
        if ~isfield(data,name)
            fprintf('   >>> Variable %s not found in existing data\n',name)
            data.(name)=vals;
            for j=1:numel(comadex)+1
                data.([name '_nmlindex' num2str(j)])=index(j);                
            end
            data.([name '_maxorder'])=numel(comadex)+1;
            fprintf('   >>> Variable %s has %d dimensions\n',name,numel(comadex)+1)
        else
            fprintf('   >>> Variable %s found in existing data\n',name)
            data.(name)=[data.(name) vals];
            for j=1:numel(comadex)+1
                data.([name '_nmlindex' num2str(j)])=[data.([name '_nmlindex' num2str(j)]) index(j)];
            end
        end
    end
end
if DEBUG 
    fprintf('------------------------------------------\n')
    fprintf(' Done parsing, now reformulating\n')
    fprintf('------------------------------------------\n')
end

% Now we reformulate the arrays
data=rmfield(data,'temp');
names=fieldnames(data);
for i=1:numel(names)
    if isfield(data,[names{i} '_nmlindex'])  % 1D array with some indices
        if DEBUG
            fprintf('Reformulating variable %s\n',names{i})
        end
        for j=1:numel(data.(names{i}))
            temp(data.([names{i} '_nmlindex'])(j))=data.(names{i})(j);
        end
        data.(names{i})=temp;
        data=rmfield(data,[names{i} '_nmlindex']);
    elseif isfield(data,[names{i} '_nmlindex1']) % multi-D array with indices
        if DEBUG
            fprintf('Reformulating variable %s\n',names{i})
        end
        if DEBUG
            fprintf('   >>> Found _nmlindex1\n')
        end
        test=0;
        for j=1:data.([names{i} '_maxorder'])
            if isfield(data,[names{i} '_nmlindex' num2str(j)])
                test=1;
            else
                test=0;
            end
        end
        if test
            minmax=zeros(2,data.([names{i} '_maxorder']));
            for j=1:data.([names{i} '_maxorder'])
                minmax(1,j)=min(data.([names{i} '_nmlindex' num2str(j)]));
                minmax(2,j)=max(data.([names{i} '_nmlindex' num2str(j)]));
            end
            arraysize=minmax(2,1)-minmax(1,1)+1;
            for j=2:data.([names{i} '_maxorder'])
                arraysize=[arraysize minmax(2,j)-minmax(1,j)+1];
            end
            if DEBUG
                fmt = repmat('%d ',[1,size(arraysize,2)]);
                fprintf(['   >>>Determined array size to be ',fmt,'\n'],arraysize)
            end
            temp=zeros(arraysize);
            % Now we redo the _nmlindex arrays to get the proper ordering
            for j=1:data.([names{i} '_maxorder'])
                if minmax(1,j)< 1
                    data.([names{i} '_nmlindex' num2str(j)])=...
                        data.([names{i} '_nmlindex' num2str(j)])-minmax(1,j)+1;
                end
                if minmax(2,j)< 1
                    data.([names{i} '_nmlindex' num2str(j)])=...
                        data.([names{i} '_nmlindex' num2str(j)])-minmax(2,j)+1;
                end                
            end
            % Create the new array by creating a string to index multiple
            % indexes
            for j=1:numel(data.(names{i}))
                % Create a executable string
                exestring=['temp('];
                for k=1:data.([names{i} '_maxorder'])                    
                    temp_ind = j;
                    if temp_ind > numel(data.([names{i} '_nmlindex' num2str(k)]))                        
                        if k == 1
                            data.([names{i} '_nmlindex' num2str(k)])(temp_ind) = data.([names{i} '_nmlindex' num2str(k)])(temp_ind-1) + 1;
                        else
                            temp_ind = 1;
                        end
                    end
                    exestring=[exestring 'data.' names{i} '_nmlindex' num2str(k)...
                        '(' num2str(temp_ind) '),'];
                end
                exestring=[exestring(1:numel(exestring)-1) ')=data.' names{i}...
                    '(' num2str(j) ');'];
                exestring
                eval(exestring);
            end
            data.(names{i})=temp;
            % Now cleanup fields
            for k=1:data.([names{i} '_maxorder'])
                data=rmfield(data,[names{i} '_nmlindex' num2str(k)]);
            end
            data=rmfield(data,[names{i} '_maxorder']);
        end
    end
end
return
end

function output=handlestars(input)
% This function handles the multiple field references made in an input
% namelist such as 128*5.5.
output=input;
stardex=strfind(output,'*');
while ~isempty(stardex)
        spdex1=strfind(output(1:stardex(1)),' ');
        spdex1=spdex1(numel(spdex1));
        numvals=sscanf(output(spdex1:stardex(1)-1),'%g');
        spdex2=strfind(output(stardex(1):numel(output)),' ');
        front=output(1:spdex1);
        if isempty(spdex2)
            spdex2=numel(output);
            val=sscanf(output(stardex(1)+1:numel(output)),'%g');
            back='';
        else
            spdex2=spdex2(1)+stardex(1);
            val=sscanf(output(stardex(1)+1:spdex2-1),'%g');
            back=output(spdex2:numel(output));
        end
        output=front;
        for i=1:numvals
            output=[output ' ' num2str(val)];
        end
        output=[output ' ' back];
        stardex=strfind(output,'*');
end
return
end