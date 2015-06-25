function write_vmec_coils_file(fname,cstruct)

fid = fopen(fname,'w');
fprintf(fid,'periods %d\n',cstruct.num_periods);
fprintf(fid,'%s\n',char(cstruct.fil_str));
fprintf(fid,'%s\n',char(cstruct.mir_str));

num_coils = length(cstruct.cname);
for i = 1:num_coils
    npts = cstruct.numelc(i);

    for j = 1:npts-1
        fprintf(fid,'%e %e %e %e\n',cstruct.coilpos{i}(j,1:3),cstruct.coilcur{i}(j));
    end
    fprintf(fid,'%e %e %e %e %d %s\n',cstruct.coilpos{i}(j,1:3),cstruct.coilcur{i}(npts),cstruct.cnum(i),char(cstruct.cname{i}));
end

fprintf(fid,'end\n');

fclose(fid);
