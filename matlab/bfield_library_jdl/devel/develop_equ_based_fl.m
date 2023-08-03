% gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\KSTAR\K-DEMO_files_mhchoi\'

clearvars;

gfile_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Haskey\SOLPS\high_n_449\baserun\gfile';
equ_name = 'C:\Users\jjl\Dropbox (ORNL)\DIII-D\Haskey\SOLPS\high_n_449\baserun\g179477.02280.equ';




R = 2.1;
Z = -.1;

Rstart = R;
Zstart = Z;
phistart = 0;
dl = 1e-2;
nsteps = 100;

figure; hold on; box on; grid on; set(gcf,'color','w');


g = readg_g3d(gfile_name);
Bout = bfield_geq_bicub(g,R,Z);
bfield.type = 'gfile';
bfield.nsym = 1;
bfield.g = g;
s = follow_fieldlines_rzphi_dl(bfield,Rstart,Zstart,phistart,dl,nsteps);
fprintf('R = %e, Z = %e, B = [%e, %e, %e]\n',R,Z,Bout.br,Bout.bz,Bout.bphi);
plot(s.r,s.z,'LineWidth',2)


clear bfield s
Equ = read_equ_file(equ_name);
Bout = bfield_equ_bicub(Equ,R,Z);
bfield.type = 'equ';
bfield.nsym = 1;
bfield.Equ = Equ;
fprintf('R = %e, Z = %e, B = [%e, %e, %e]\n',R,Z,Bout.br,Bout.bz,Bout.bphi);
s = follow_fieldlines_rzphi_dl(bfield,Rstart,Zstart,phistart,dl,nsteps);
plot(s.r,s.z,'--','LineWidth',2)
