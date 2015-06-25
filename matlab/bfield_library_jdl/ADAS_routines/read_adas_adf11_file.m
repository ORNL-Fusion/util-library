function [te,ne,cc] = read_adas_adf11_file(fname)
% te = eV, ne = cm^-3, cc = cc(iz,ne,te)
% coeff = zeros(iz1max,idmaxd,itmaxd);
% clearvars; 

% type = 'scd';

% fname = ['C:\Work\ADAS\ADF11\','scd96_h.dat'];
% fname = ['C:\Work\ADAS\adf11_all\scd96\','scd96_n.dat'];

fid = fopen(fname,'r');
if fid == -1    
    error(['Did not find output file: ',fname])
end

% First line
nums = fscanf(fid,'%i %i %i %i %i'); izmax = nums(1); idmaxd = nums(2); itmaxd = nums(3); iz1min = nums(4); iz1max = nums(5);
strs = fscanf(fid,'%s',1); species_name = strs(2:end);
strs = fgetl(fid); note = strs(2:end);

fprintf('Reading file: %s\n',fname);
fprintf('Species: %s\n',species_name);
fprintf('Nuclear charge %i\n',izmax);
fprintf('Number of densities, temperatures = [%i,%i]\n',idmaxd,itmaxd);
fprintf('Lowest, highest charge = [%i,%i]\n',iz1min,iz1max);

% Read ne, Te
junk = fgetl(fid);
ddensd = fscanf(fid,'%f\n',idmaxd);  %log10(ne) (cm-3)
dtevd = fscanf(fid,'%f\n',itmaxd);   %log10(Te) (eV)

% Loop over charge state and read collisional radiative coefficients
% coeff = zeros(iz1max,idmaxd,itmaxd);
for iz = 1:iz1max
    junk = fgetl(fid);  % could mine this a bit
    coeff_tmp = fscanf(fid,'%f\n',itmaxd*idmaxd);
    coeff(iz,:,:) = reshape(coeff_tmp,idmaxd,itmaxd);  % Should be log10 of coefficient    
end
% 
% figure; hold on; box on;
% for iz = 1:iz1max
%     plot(dtevd,squeeze(coeff(iz,:,:)))
% end

ne = 10.^ddensd;
te = 10.^dtevd;
cc = 10.^coeff;
% figure; hold on; box on;
% for iz = 1:iz1max
% %     plot(te,squeeze(cc(iz,:,:)))
%     plot(te,squeeze(cc(iz,1,:)))
% end
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('T_e (eV)','fontsize',12)
% ylabel('Coeff (cm^3/s)','fontsize',12)
% set(gca,'fontsize',12)


