function ipec = open_ipec_field(run_path)

% clearvars;

% run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\high\gpec\';
% run_path = 'C:\Work\DIII-D\164723\g164723.03059_d3d_kinetic\low\gpec\';

fname_eq = fullfile(run_path,'ipec_eqbrzphi_n3.out');
fname_p  = fullfile(run_path,'ipec_pbrzphi_n3.out');
fname_c  = fullfile(run_path,'ipec_cbrzphi_n3.out');

ipec.eq = read_ipec_field_file(fname_eq);
ipec.pert = read_ipec_field_file(fname_p);
ipec.vac = read_ipec_field_file(fname_c);



% figure; hold on; box on;
% plot(ipec.eq.r,ipec.eq.z,'k.')
% plot(ipec.pert.r,ipec.pert.z,'ro')
% plot(ipec.vac.r,ipec.vac.z,'gx')