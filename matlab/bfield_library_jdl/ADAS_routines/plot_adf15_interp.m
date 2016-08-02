% function plot_adf15_interp
clearvars;

% wl_nm = 656;
% fname = 'C:\Work\ADAS\pec12#h_pju#h0.dat';

wl_nm = 514;
fname = 'C:\Work\ADAS\pec96#c_vsu#c1.dat';

adf15 = read_adas_adf15_file(fname);
[i_ex,i_rc] = find_adf15_transition(adf15,wl_nm,1);

ne_test_m3 = 1e20;
te_test_eV = logspace(log10(1),log10(2000),50);

sigma_v_ex = interp_adas_rate_coefficient(te_test_eV,ne_test_m3/1e6,adf15.te{i_ex},adf15.dens{i_ex},adf15.pec{i_ex});
sigma_v_rc = interp_adas_rate_coefficient(te_test_eV,ne_test_m3/1e6,adf15.te{i_rc},adf15.dens{i_rc},adf15.pec{i_rc});

plot(te_test_eV,sigma_v_ex,'r.')
plot(te_test_eV,sigma_v_rc,'b.')