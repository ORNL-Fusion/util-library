clearvars;

gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Qprl experiment\power13mw\g174308.03500_159';
g = readg_g3d(gfile_name);


Rtest = linspace(1.8,3.4,100);
Ztest = linspace(-1,1,100);



bb = bfield_geq_bicub_inv(g,Rtest,Ztest);
% fprintf('R, Z, Br, Bphi, Bz = [%8.6f, %8.6f, %8.6f, %8.6f, %8.6f]\n',Rtest,Ztest,bb.br,bb.bphi,bb.bz)

