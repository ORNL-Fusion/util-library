clearvars;

% url_base = 'http://open.adas.ac.uk/detail/adf11/';
url_base = 'http://open.adas.ac.uk/download/adf11/';


types = {'acd','plt','scd','ccd','prb','prc'};
species = {'c','h','he','li','be','ar','o','ne','n','w','fe','xe'};
years = {'89','93','96','12'};

outname_base = 'C:\Work\ADAS\adf11_all';

for k = 1:length(species)
    for j = 1:length(years)
        for i = 1:length(types)
            url_name = strcat(url_base,types{i},years{j},'/',types{i},years{j},'_',species{k},'.dat')
            outdir = fullfile(outname_base,strcat(types{i},years{j}));
            mkdir(outdir);
            outfile = fullfile(outdir,strcat(types{i},years{j},'_',species{k},'.dat'));
            test = websave(outfile,url_name);
            [a,b,c] = fileparts(test);
            if ~strcmp(c,'.dat')
                delete(test)
            end
            test2 = dir(outdir);
            if length(test2) == 2  % empty
                rmdir(outdir)
            end
%             urlwrite(url_name,outfile);
        end
    end
end

