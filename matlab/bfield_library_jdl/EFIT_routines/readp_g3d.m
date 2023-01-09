function p = readp_g3d(filename)

% filename = 'C:\Work\DIII-D\160884\efits\p160884.03014_251';

fid = fopen(filename,'r');

counter = 1; counter_max = 1000;

%%
while ~feof(fid)
    line = fgetl(fid);

    %% Handle "N Z A ION SPECIES" line
    if contains(line,'ION SPECIES','IgnoreCase',true)
        p.nions = sscanf(line,'%d');
        p.ion_data_NZA = reshape(fscanf(fid,'%f %f %f\n',p.nions*3),3,p.nions).';
        continue;
    end

    %% Otherwise parse based on string for profile data
    % Get npts
    p.npts = sscanf(line,'%d');
    if counter > 1
        if p.npts ~= npts_last
            error('Different value for npts, have to account for this')
        end
    end
    npts_last = p.npts;
    % Reshape data
    temp = reshape(fscanf(fid,'%f %f %f\n',p.npts*3),3,p.npts).';
    this = line(strfind(line,'psinorm')+8:strfind(line,"(")-1);
    % Check string and assign
    p.profiles.(this).psin = temp(:,1);
    p.profiles.(this).data = temp(:,2);
    p.profiles.(this).dpsin = temp(:,3);
    p.profiles.(this).units = line(strfind(line,"(")+1:strfind(line,")")-1);
end
fclose(fid);

counter = counter + 1;
if counter > counter_max
    error('Hit counter_max in readp_g3d')
end


%% Check if psi is uniform
names = fieldnames(p.profiles);
p1 = p.profiles.(names{1}).psin;
uniform = 1;
for i = 2:length(names)
    if ~all(p1 == p.profiles.(names{i}).psin)
        uniform = 0;
    end
end
if ~uniform
    fprintf('psiN is not uniform across profiles!\n');
else
    fprintf('psiN is uniform across profiles, consolidating\n');
    p.profiles.psin = p1;
    for i = 1:length(names)
        p.profiles.(names{i}) = rmfield(p.profiles.(names{i}),'psin');
    end

end





