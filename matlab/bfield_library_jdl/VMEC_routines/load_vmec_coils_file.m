function [coil,current,cnum,cname,numelc,cstruct]=load_vmec_coils_file(fname)
% fname = 'coils.w7x';
max_count = 1e6;

% Get number of lines in file
fid = fopen(fname,'r');
s = fgetl(fid);
s = fgetl(fid);
s = fgetl(fid);
num_lines = 0;
while ~strcmp(s,'end') 
    s = fgetl(fid);
    num_lines = num_lines + 1;
    if num_lines > max_count 
        error('max count exceeded in coils file')
    end
end
num_lines = num_lines -1;
fclose(fid);

% Parse file
fid = fopen(fname,'r');
data = fgetl(fid); data2 = textscan(data,'%s %d'); num_periods = data2{2}(1);
data = fgetl(fid); data2 = textscan(data,'%s'); fil_str = [char(data2{1}{1}),' ',char(data2{1}{2})];
data = fgetl(fid); data2 = textscan(data,'%s'); mir_str = [char(data2{1}{1}),' ',char(data2{1}{2})];

current = zeros(num_lines,1);
coil = zeros(num_lines,3);
ccount = 0;
numel = 0;

% for structure
icount_intern = 1;
icount_coil = 1;

for i = 1:num_lines
    data = fgetl(fid);
    [~,num_entries]=sscanf(data,'%s');
    numel = numel + 1;
    if num_entries == 4
        data2 = textscan(data,'%f %f %f %f');      
    else
        data2 = textscan(data,'%f %f %f %f %d %s');
        ccount = ccount + 1;
        cnum(ccount) = data2{5}(1);
        cname{ccount} = data2{6};
        numelc(ccount) = numel;
        numel = 0;
    end
    coil(i,1:3) = [data2{1}(1) data2{2}(1) data2{3}(1)];
    current(i) = data2{4}(1);
    
    % separated by coil
    cstruct.coilpos{icount_coil}(icount_intern,1:3) = coil(i,1:3);
    cstruct.coilcur{icount_coil}(icount_intern)     = current(i);
    if numel == 0
        icount_coil = icount_coil + 1;
        icount_intern = 0;
    end
    icount_intern = icount_intern + 1;
end


fclose(fid);

if nargout > 5
    cstruct.coil = coil;
    cstruct.current = current;
    cstruct.cnum = cnum;
    cstruct.cname = cname;
    cstruct.numelc = numelc;
    cstruct.num_periods = num_periods;
    cstruct.fil_str = fil_str;
    cstruct.mir_str = mir_str;
end

end
