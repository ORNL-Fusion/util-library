function [mfp,sigma_v] = calc_mfp(te,ne,RateCoeff,indReac,vp)
% Inputs:
%   Te in eV
%   ne in m^-3
%   RateCoeff from read_adas_adf11_file (or equivalent) cm^3/s   (Wcm^3 for plt)
%   indReac: value for last index in RateCoeff.coeff (if scd, ind = charge + 1)
%   vp: test particle velocity (m/s)
% Outputs:
%   mfp in m


ind_oob = find(ne > max(RateCoeff.ne)*1e6);
if ~isempty(ind_oob)
    warning('Some densities are out of bounds in upper direction! Max value used!!')
    ne(ind_oob) = max(RateCoeff.ne)*1e6;
end

sigma_v = 1e-6*10.^(interp2(RateCoeff.te_log10,RateCoeff.ne_log10+6,squeeze(RateCoeff.coeff_log10(:,:,indReac)),log10(te),log10(ne),'linear',NaN));
mfp = vp./(ne.*sigma_v);
