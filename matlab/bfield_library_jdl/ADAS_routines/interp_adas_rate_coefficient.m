function sigma_v = interp_adas_rate_coefficient(Te,ne,te_adas,ne_adas,coeff_adas)
% Te in eV, ne in cm^-3
% coefficient cm^3/s  

log_interp = 1;
my_method = 'spline';
% my_method = 'linear';

ind_oob = find(ne > max(ne_adas));
if ~isempty(ind_oob)
    warning('Some densities are out of bounds!')
    ne(ind_oob) = max(ne_adas);
end


extrap_val = NaN;
if log_interp == 1
    sigma_v = exp(interp2(log(te_adas),log(ne_adas),log(coeff_adas),log(Te),log(ne),my_method,extrap_val));    
else
    sigma_v = interp2(te_adas,ne_adas,coeff_adas,Te,ne,my_method,extrap_val);
end

