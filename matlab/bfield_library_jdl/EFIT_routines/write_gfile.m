% function write_gfile(g,fname)


gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\cd-46864-400ms.geqdsk';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\ed-47079-400ms.geqdsk';
% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\sxd-46860-500ms.geqdsk';
g = readg_g3d_simple(gfile_name);


fname = 'C:\Users\jjl\Dropbox (ORNL)\MAST\gfiles_from_jack\cd-46864-400ms_formatted.geqdsk';
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

%%
fprintf(fid,f2020,g.xdim,g.zdim,g.rzero,g.rgrid1,g.zmid);
fprintf(fid,f2020,g.rmaxis,g.zmaxis,g.ssimag,g.ssibry,g.bcentr);
fprintf(fid,f2020,g.cpasma,g.ssimag,xdum,g.rmaxis,xdum);
fprintf(fid,f2020,g.zmaxis,xdum,g.ssibry,xdum,xdum);
fprintf(fid,f2020,g.fpol);
fprintf(fid,f2020,g.pres);
fprintf(fid,f2020,g.ffprim);
fprintf(fid,f2020,g.pprime);
fprintf(fid,f2020,g.psirz);
if mod(numel(g.psirz),5) ~= 0 
    fprintf(fid,'\n');
end
fprintf(fid,f2020,g.qpsi);
fprintf(fid,f2022,g.nbdry,g.limitr);
fprintf(fid,f2020,g.bdry);
if mod(numel(g.bdry),5) ~= 0 
    fprintf(fid,'\n');
end
fprintf(fid,f2020,g.lim);
fprintf(fid,'\n');

fclose(fid);


