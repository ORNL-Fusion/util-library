% RESTORE_IDL_DEMO
%
% finds the IDL save files included in the restore_idl distribution and
% restores each of them in turn

me = evalc('which restore_idl_demo');
[path, file, type]=fileparts(me);
savefiles=dir(fullfile(path,'idl_save','*.sav'));
for i=1:numel(savefiles)
    fname=savefiles(i).name;
    disp(['press any key to restore ' fname]);
    pause;
    % restore_idl places the variables into a structure with field names
    % that are the IDL variable names, field values are the IDL variables
    % themselves.
    outargs=restore_idl(fullfile(path,'idl_save',fname),'verbose');
    outargs;
    % show what we found
    disp('Contents of outargs:');
    fields=fieldnames(outargs)
    for n=1:numel(fields)
        % use eval to create a variable from the structure field
        eval([fields{n} '=outargs.' fields{n} ';']);
        % if number of elements is small, show them
        if numel(outargs.(fields{n}))<=10,
            eval((fields{n}))
        else
            eval(['whos ' fields{n}])
        end
    end
end