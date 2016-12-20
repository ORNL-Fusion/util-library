function query_vmec_coils_file(coils_file)
if ~isstruct(coils_file)
    coil = load_vmec_coils_file(coils_file);
else
    coil = coils_file;
end

fprintf('\n--------------------------------------------\n')
fprintf('All coil information\n')
fprintf('--------------------------------------------\n')
for i = 1:length(coil.cname)
    fprintf('Coil number %3d, index %2d, factor %8.2f, %s\n',i,coil.cnum(i),coil.numwind(i),char(coil.cname{i}))
end

[inds,ia] = unique(coil.cnum);
fprintf('\n--------------------------------------------\n')
fprintf('Unique coil information\n')
fprintf('--------------------------------------------\n')
for i = 1:length(ia)
    fprintf('Coil index %2d, factor %8.2f, %s\n',coil.cnum(ia(i)),coil.numwind(ia(i)),char(coil.cname{ia(i)}))
end