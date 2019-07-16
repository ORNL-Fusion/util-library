% function TEST_g3d_eval
clearvars;

gfile_name = 'C:/Work/DIII-D/148712/g148712.04101';

g = readg_g3d(gfile_name);

N = 10000000;
a = 1; b = 2.4;
R = sort(a + (b-a).*rand(N,1));

a = -1.2; b = 1.2;
Z = a + (b-a).*rand(N,1);


% profile on;
tic; 
Binv = bfield_geq_bicub_inv(g,R,Z);
toc
% profile report

tic;
% profile on;
Borig = bfield_geq_bicub(g,R,Z);
% profile report
toc
% 
% figure; hold on; 
% plot(R,Borig.br)
% plot(R,Binv.br,'--')
% 
% figure; hold on; 
% plot(R,B.br-B2.br)