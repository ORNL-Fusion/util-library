function fname = find_IR_file(shot)
% Search data path for file names including string of shot number
data_path = 'C:\Work\MPEX\Camera\';
files = dir(data_path);
icount = 0;
for i = 3:length(files)
    itmp = strfind(files(i).name,num2str(shot));
    if ~isempty(itmp)
        ishot = i;
        icount = icount + 1;   
    end
end
if icount == 0
    error(['Could not find camera data for shot ',num2str(shot)])
end
if icount > 1
    error(['Found multiple files matching shot: ',num2str(shot)])
end
fname = strcat(data_path,files(ishot).name);
fprintf('Using file %s \n',fname)