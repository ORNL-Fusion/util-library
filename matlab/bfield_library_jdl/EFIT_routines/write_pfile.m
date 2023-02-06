function write_pfile(p,fileName)

fid = fopen(fileName,'w');

ProfNames = fieldnames(p.profiles);

if isfield(p.profiles,'psin')
    uniformPsiN = 1;
    psiN = p.profiles.psin;
else
    uniformPsiN = 0;
end

%% Write profiles
for i = 1:length(ProfNames)  - uniformPsiN 
    thisName = ProfNames{i};
    this = p.profiles.(thisName);
    fprintf(fid,'%d %s %s %s\n',length(this.data),'psinorm',strcat(thisName,'(',this.units,')'),strcat('d',thisName,'/dpsiN'));
    if ~uniformPsiN
        psiN = this.psiN;
    end
    fprintf(fid,'%9.6f %10.6f %11.6f\n',[psiN,this.data,this.dpsin]');
end
%% Write "N Z A ION SPECIES" line
if isfield(p,'ion_data_NZA')
    fprintf(fid,'%d N Z A of ION SPECIES\n',p.nions);
    fprintf(fid,'%9.6f %9.6f %9.6f\n',p.ion_data_NZA');
end

fclose(fid);