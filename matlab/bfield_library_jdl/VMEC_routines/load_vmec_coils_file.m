function cstruct=load_vmec_coils_file(fname)
% old output:[coil,current,cnum,cname,numelc,cstruct]
% fname = 'coils.w7x';

% DISPLAY_COIL_NAMES = 1;

num_lines = num_lines_file(fname);

% Parse file
fprintf('Reading coils file %s\n',fname)
fid = fopen(fname,'r');
data = fgetl(fid);
if ~strcmp(data(1:7),'periods')
max_count = 100;  % Read MGRID_MLI namelist if it exists
for i0 = 1:max_count
    if i0 > 1
        data = fgetl(fid); 
    end
    if strcmp(data,'** coils_dot_starts_below **')
        break;
    end
end
data = fgetl(fid);
end
data2 = textscan(data,'%s %f'); num_periods = data2{2}(1);
fprintf('Number of periods %f\n',num_periods);
data = fgetl(fid); data2 = textscan(data,'%s'); fil_str = [char(data2{1}{1}),' ',char(data2{1}{2})];
data = fgetl(fid); data2 = textscan(data,'%s'); mir_str = [char(data2{1}{1}),' ',char(data2{1}{2})];

% current = zeros(num_lines,1);
% coil = zeros(num_lines,3);
ccount = 0;
numel = 0;

% for structure
icount_intern = 1;
icount_coil = 1;

for i = 1:num_lines-i0
    data = fgetl(fid);
    if strcmp(data,'end')
        break;
    end
    [~,num_entries]=sscanf(data,'%s');
    numel = numel + 1;
    if num_entries == 4
        data2 = textscan(data,'%f %f %f %f');      
        numwind(ccount+1) = data2{4}(1);
    else        
        data2 = textscan(data,'%f %f %f %f %f %s');
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


    cstruct.coil = coil;
    cstruct.current = current;
    cstruct.cnum = cnum;
    cstruct.cname = cname;
    cstruct.numelc = numelc;
    cstruct.num_periods = num_periods;
    cstruct.fil_str = fil_str;
    cstruct.mir_str = mir_str;
    cstruct.numwind = numwind;


end
