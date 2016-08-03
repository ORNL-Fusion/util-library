function field = read_ipec_field_file(fname)
% field.r(ir,iz) -- indices start from lower left corner

% run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\high\gpec\';
% run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\low\gpec\';
% fname = fullfile(run_path,'ipec_eqbrzphi_n3.out');
% fname = fullfile(run_path,'ipec_cbrzphi_n3.out');


[base,fname_part,ext] = fileparts(fname);
matfile = fullfile(base,[fname_part,ext,'.mat']);
if exist(matfile,'file')
    fprintf('Reading ipec field data from .mat file\n')
    load(matfile);
    return;
else
    fprintf('Reading ipec data from ASCII file\n')
end

fid = fopen(fname);
dat = fgetl(fid);
fprintf('%s\n',dat)
if ~isempty(strfind(dat,'_EQBRZPHI'))
    field_size = 3;
elseif ~isempty(strfind(dat,'_PBRZPHI'))
    field_size = 6;
elseif ~isempty(strfind(dat,'_CBRZPHI'))
    field_size = 6;
elseif ~isempty(strfind(dat,'_BRZPHI'))
    field_size = 6;    
else
    error('Did not recognize field type: %s\n',dat)
end

dat = fgetl(fid); % version
dat = fgetl(fid); % blank
dat = fgetl(fid);  % n = 
dat2 = sscanf(dat,'%s %s %d');
n = dat2(3);
fprintf('n = %d\n',n)

dat = fgetl(fid);  % nr, nz
dat2 = sscanf(dat,'%s %s %d %s %s %d');
nr = dat2(4); nz = dat2(8);
fprintf('nr = %d, nz = %d\n',nr,nz)

dat = fgetl(fid); % blank
dat = fgetl(fid);  % label

if field_size == 3
    data = fscanf(fid,'%d %f %f %f %f %f\n',[3+field_size,nr*nz]);
else
    data = fscanf(fid,'%d %f %f %f %f %f %f %f %f\n',[3+field_size,nr*nz]);
end

fclose(fid);


field.n = n;
field.nr = nr;
field.nz = nz;
field.npts = nr*nz;
field.r = reshape(data(2,:),[nr,nz]).';
field.z = reshape(data(3,:),[nr,nz]).';
if field_size == 3    
    field.br   = reshape(data(4,:),[nr,nz]).';
    field.bz   = reshape(data(5,:),[nr,nz]).';
    field.bphi = reshape(data(6,:),[nr,nz]).';
else
    field.rbr = reshape(data(4,:),[nr,nz]).';
    field.ibr = reshape(data(5,:),[nr,nz]).';
    field.rbz = reshape(data(6,:),[nr,nz]).';
    field.ibz = reshape(data(7,:),[nr,nz]).';
    field.rbphi = reshape(data(8,:),[nr,nz]).';
    field.ibphi = reshape(data(9,:),[nr,nz]).';
end

save(matfile,'field')
