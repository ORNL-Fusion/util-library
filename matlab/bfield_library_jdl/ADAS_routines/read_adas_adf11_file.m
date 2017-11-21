function data = read_adas_adf11_file(fname)
% reads scd, acd, plt, prb, ccd, etc
%
% Note: Old call:
%   [te,ne,cc] = read_adas_adf11_file(fname)
%   te = eV, ne = cm^-3, cc = cc(iz,ne,te)
%   old cc was (iz,ne,te), new rate is (ne,te,iz)
%   

fid = fopen(fname,'r');
if fid == -1    
    error('Did not find output file: %s\n',fname)
end

% First line
nums = fscanf(fid,'%d %d %d %d %d'); izmax = nums(1); idmaxd = nums(2); itmaxd = nums(3); iz1min = nums(4); iz1max = nums(5);
strs = fscanf(fid,'%s',1); species_name = strs(2:end);
strs = fgetl(fid); note = strs(2:end);

if izmax ~= iz1max
    error('Have to update this')
end

fprintf('Reading %s file: %s\n',species_name,fname);
% fprintf('    Species, nuclear charge: %s, %d\n',species_name,izmax);
% fprintf('    Nuclear charge %d\n',izmax);
% fprintf('    Number of densities, temperatures = [%d,%d]\n',idmaxd,itmaxd);
% fprintf('    Lowest, highest charge = [%d,%d]\n',iz1min,iz1max);

% Read ne, Te
junk = fgetl(fid);
ddensd = fscanf(fid,'%f\n',idmaxd);  %log10(ne) (cm-3)
dtevd = fscanf(fid,'%f\n',itmaxd);   %log10(Te) (eV)

% Loop over charge state and read collisional radiative coefficients
% coeff = zeros(idmaxd,itmaxd,iz1max);
for iz = 1:iz1max
    junk = fgetl(fid);  % could mine this a bit
    coeff_tmp = fscanf(fid,'%f\n',itmaxd*idmaxd);
    coeff(:,:,iz) = reshape(coeff_tmp,idmaxd,itmaxd);  % Should be log10 of coefficient    
end
fclose(fid);
% 
% figure; hold on; box on;
% for iz = 1:iz1max
%     plot(dtevd,squeeze(coeff(iz,:,:)))
% end

for iz = 1:iz1max
    [tt,nn] = ndgrid(dtevd,ddensd);
    data.interp_log10(iz).interpolant = griddedInterpolant(tt,nn,coeff(:,:,iz).','spline');
end
    



data.ne_log10    = ddensd;
data.te_log10    = dtevd;
data.coeff_log10 = coeff;
data.ne          = 10.^ddensd;
data.te          = 10.^dtevd;
data.coeff       = 10.^coeff;
data.charge      = iz1max;
% data.interpolant_log10 = interpolant_log10;

% figure; hold on; box on;
% for iz = 1:iz1max
% %     plot(te,squeeze(coeff(:,:,iz)))
%     plot(te,squeeze(coeff(1,:,iz)))
% end
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlabel('T_e (eV)','fontsize',12)
% ylabel('Coeff (cm^3/s)','fontsize',12)
% set(gca,'fontsize',12)


