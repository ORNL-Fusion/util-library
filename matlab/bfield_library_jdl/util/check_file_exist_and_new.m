function [exists,new] = check_file_exist_and_new(fname1,fname2,iwarn,purge)
    % Check if fname1 exists (output as logical "exists")
    % If it exists, check if fname2 exists. If it does, output "new" is set to 1 if
    % fname1 is newer or the same timestamp as fname2. If fname2 does not exist new = 1;
    % if purge == 1 then fname2 is deleted
    % Example call: [exists,new] = check_file_exist_and_new("mydata.raw","mydata.mat",1,0)
    
    new = 1;
    exists = 1;
    if exist(fname1,'file') ~= 2
        if iwarn
            fprintf('Did not find output file: %s\n',fname1)
        end
        exists = 0;
        return;
    end
    if exist(fname2,'file') == 2
        if purge == 1
            delete(fname2);
            return;
        end
        finfo1 = dir(fname1);
        finfo2 = dir(fname2);
        if datenum(finfo2.date) > datenum(finfo1.date)
            new = 0;
        end
    end
end