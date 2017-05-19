clearvars;
run_path = 'C:\Work\FLARE\2dts_scan_as_4000ms\';
nz = 3;

d_vac_raw_2 = dlmread([run_path,'debug/VacuumBoundary_2_8.raw']);
C_cutL = dlmread([run_path,'C_cutL.plt']);
C_cutR = dlmread([run_path,'C_cutR.plt']);
guiding = dlmread([run_path,'wall_align.dat']);
rpath_2 = dlmread([run_path,'rpath_2.plt']);

figure; hold on; box on;
plot(d_vac_raw_2(:,1),d_vac_raw_2(:,2),'r','linewidth',2)
plot(C_cutL(:,1),C_cutL(:,2),'k-.','linewidth',2)
plot(C_cutR(:,1),C_cutR(:,2),'k--','linewidth',2)
plot(guiding(:,1),guiding(:,2),'k-','linewidth',2)
plot(rpath_2(:,1),rpath_2(:,2),'b-','linewidth',2)
legend('VacuumBoundary','C\_cutL','C\_cutR','guiding','rpath')