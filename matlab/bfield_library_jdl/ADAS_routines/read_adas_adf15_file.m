function adf15 = read_adas_adf15_file(fname,verbose)
% te = eV, ne = cm^-3, PEC(DEN,TEMP)  (cm3/s)
% wlng = ang
% clearvars; 
if nargin < 2
    verbose = 1;
end

% create structure
adf15 = [];

% fname = 'C:\Work\ADAS\pec12#h_pju#h0.dat';

fid = fopen(fname,'r');
if fid == -1    
    error(['Did not find output file: ',fname])
end

% First line
line = fgetl(fid); 
adf15.nsel = sscanf(line(1:5),'%d5'); % NSEL = number of transitions
adf15.text = line(6:end);   
if verbose
    fprintf('Reading ADAS15 file: %s\n',fname)
    fprintf('Label: %s\n',adf15.text)
    fprintf('Number of transitions = %d\n',adf15.nsel)
end

for isel = 1:adf15.nsel
    line = fgetl(fid); data = sscanf(line,'%fA %d %d');
    adf15.wlng(isel)  = data(1);   % wavelength of transition (Ang)
    adf15.ndens(isel) = data(2);   % number of densities
    adf15.nte(isel)   = data(3);   % number of temperatures
    sis = find(line == '/');
    if length(sis) ~= 4
        error('Could not parse adf15 file (slashes)')
    end
    for i = 1:length(sis)        
        if i < length(sis)
            i2 = sis(i+1)-1;
        else
            i2 = length(line);
        end
        switch i            
            case 1
                check = 'FILMEM';    % Source specific ion excitation file
                fmt = '%s';                    
            case 2
                check = 'TYPE';      % Type of photon emissivity (excit, recomb, cx)
                fmt = '%s';
            case 3
                check = 'INDM';      % Associated metastable index in source file
                fmt = '%d';
            case 4
                check = 'ISEL';      % Transition index
                fmt = '%d';                
        end
        if ~strcmp(line(sis(i)+1:sis(i)+length(check)),check)
            error(['Could not parse adf15 file:',check])
        end
        i1 = sis(i)+length(check)+3;
        switch i            
            case 1
                adf15.filemem{isel} = sscanf(line(i1:i2),fmt);
            case 2
                adf15.type{isel} = sscanf(line(i1:i2),fmt);
            case 3
                adf15.indm(isel) = sscanf(line(i1:i2),fmt);
            case 4
                adf15.isel(isel) = sscanf(line(i1:i2),fmt);
        end        
        
    end
    adf15.dens{isel} = fscanf(fid,'%f\n',adf15.ndens(isel));   % ne (cm-3)
    adf15.te{isel} = fscanf(fid,'%f\n',adf15.nte(isel));       % Te (eV)
    adf15.pec{isel} = fscanf(fid,'%f\n',[adf15.nte(isel),adf15.ndens(isel)]).';  % PEC(DEN,TEMP)  (cm3/s)     
end

% figure; hold on; box on;
% for i = 1:adf15.nsel
%     plot(adf15.te{isel},adf15.pec{isel})
%     set(gca,'xscale','log')
%     set(gca,'yscale','log')
% end


