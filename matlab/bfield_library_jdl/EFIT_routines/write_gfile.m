function write_gfile(g,fname)
% write_gfile(g,fname)


% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\cd-46864-400ms.geqdsk';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\ed-47079-400ms.geqdsk';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\sxd-46860-500ms.geqdsk';
% g = readg_g3d_simple(gfile_name);


% fname = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\cd-46864-400ms_formatted.geqdsk';
% fname = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\ed-47079-400ms_formatted.geqdsk';
% fname = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\sxd-46860-500ms_formatted.geqdsk';

fid = fopen(fname,'w');

xdum = 0;
f2022 = '%5d%5d\n';
f2020 = '%16.8e%16.8e%16.8e%16.8e%16.8e\n';

%% 
txt = g.line1(1:min(length(g.line1),29));
this = sprintf('%48s%4d%4d%4d',pad(txt,48,'right',' '),0,g.mw,g.mh);
% fprintf(fid,'%s\n',g.line1);
fprintf(fid,'%s\n',this);

%% Sanity checks
if size(g.lim,2) ~= g.limitr
    error('Size of lim structure and limitr do not match: %d, %d',size(g.lim,2), g.limitr)
end
if size(g.bdry,2) ~= g.nbdry
    error('Size of bdry structure and nbdry do not match: %d, %d',size(g.bdry,2), g.nbdry)
end


%%
write2020(fid,[g.xdim,g.zdim,g.rzero,g.rgrid1,g.zmid]);
write2020(fid,[g.rmaxis,g.zmaxis,g.ssimag,g.ssibry,g.bcentr]);
write2020(fid,[g.cpasma,g.ssimag,xdum,g.rmaxis,xdum]);
write2020(fid,[g.zmaxis,xdum,g.ssibry,xdum,xdum]);
write2020(fid,g.fpol)
write2020(fid,g.pres);
write2020(fid,g.ffprim);
write2020(fid,g.pprime);
write2020(fid,g.psirz);
write2020(fid,g.qpsi);
fprintf(fid,f2022,g.nbdry,g.limitr);
write2020(fid,g.bdry);
write2020(fid,g.lim);


fclose(fid);

end

function write2020(fid,this)
for i = 1:5:numel(this)
    endIdx = min(i+4, numel(this)); % Determine the end index for the current line
    fprintf(fid, '% 14.9e', this(i:endIdx));
    fprintf(fid, '\n'); % Add a newline character after each line
end
end

function write2022(fid,this)
for i = 1:5:numel(this)
    endIdx = min(i+4, numel(this)); % Determine the end index for the current line
    fprintf(fid, '%5d', this(i:endIdx));
    fprintf(fid, '\n'); % Add a newline character after each line
end
end