clearvars;

% url_base = 'http://open.adas.ac.uk/detail/adf11/';
url_base = 'https://open.adas.ac.uk/download/adf11/';
            % https://open.adas.ac.uk/download/adf11/scd89/scd89_w.dat
types = {'scd'};

% types = {'acd','plt','scd','ccd','prb','prc'};
species = {'c','h','he','li','be','ar','o','ne','n','w','fe','xe','w'};
% species = {'w'};
% years = {'89','93','96','12','22'};
% years = {'89','93','96'};

outname_base = 'C:\Users\jjl\ORNL Dropbox\Jeremy Lore\ADAS\adf11_all';

for k = 1:length(species)
    for j = 1:length(years)
        for i = 1:length(types)
            url_name = strcat(url_base,types{i},years{j},'/',types{i},years{j},'_',species{k},'.dat');
            fprintf('Trying url: %s',url_name)
            outdir = fullfile(outname_base,strcat(types{i},years{j}));
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end
            outfile = fullfile(outdir,strcat(types{i},years{j},'_',species{k},'.dat'));
            test = websave(outfile,url_name);
            [a,b,c] = fileparts(test);
            if ~strcmp(c,'.dat')
                delete(test);
            end
            if exist(test) == 2
                fprintf(' -- success\n');
            else
                fprintf(' -- fail\n');
            end
            test2 = dir(outdir);
            if length(test2) == 2  % empty
                try
                    rmdir(outdir)
                end
            end
%             urlwrite(url_name,outfile);
        end
    end
end


