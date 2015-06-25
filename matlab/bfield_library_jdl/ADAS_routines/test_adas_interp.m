% test_adas_interp
clearvars;

fname = ['C:\Work\ADAS\adf11_all\scd96\','scd96_h.dat'];  % Effective ionization coefficients (cm^-3/s) 
[te_scd,ne_scd,scd] = read_adas_adf11_file(fname);
fname = ['C:\Work\ADAS\adf11_all\acd96\','acd96_h.dat'];  % Effective recombination coefficients (cm^-3/s)
[te_acd,ne_acd,acd] = read_adas_adf11_file(fname);

Te = [1,10]
ne = [1e14,12e14]
sigma_v = interp_adas_rate_coefficient(Te,ne,te_scd,ne_scd,squeeze(scd(1,:,:)).' )