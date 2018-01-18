clearvars;

suffix = '96'; element='h';

adf11_path = 'C:\Work\ADAS\adf11_all';

fprintf('Reading scd\n');
fname = [adf11_path,'\scd',suffix,'\','scd',suffix,'_',element,'.dat'];  % Effective ionization coefficients (cm^-3/s)
scd = read_adas_adf11_file(fname);
fprintf('Reading acd\n');
fname = [adf11_path,'\acd',suffix,'\','acd',suffix,'_',element,'.dat'];  % Effective recombination coefficients (cm^-3/s)
acd = read_adas_adf11_file(fname);
fprintf('Reading plt\n');
fname = [adf11_path,'\plt',suffix,'\','plt',suffix,'_',element,'.dat'];  % Radiated power (W cm^3)
plt = read_adas_adf11_file(fname);
fprintf('Reading ccd\n');
fname = [adf11_path,'\ccd',suffix,'\','ccd',suffix,'_',element,'.dat'];  %
ccd = read_adas_adf11_file(fname);




Te = linspace(0.1,1000,1000);
ne = 1e12;

rate = scd;

for i = 1:length(Te)
    sigma_v(i) = interp_adas_rate_coefficient(Te(i),ne,rate.te,rate.ne,rate.coeff)/1e6;
end

figure; hold on;
plot(Te,sigma_v)