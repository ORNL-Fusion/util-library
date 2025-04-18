function val = eval_AMJUEL_H4_fit(H4data,ne_cm3,Te_eV)
% Evaluate H.4 AMJUEL fit
% density in cm^-3 (scalar)
% temperature in eV (can be vector)
% Rate coefficient (val) in cm^3/s
% H4data comes from function return_AMJUEL_data
%
% Warning: No sanity checks on range of ne, Te where fit is valid!!
%
% J.D. Lore
%
% Example:
% data = return_AMJUEL_data;
% Te_test = logspace(-1,2,100); % eV
% ne_test = 1e13; % cm3
% rate = eval_AMJUEL_H4_fit(data.H4.reaction_215,ne_test,Te_test);
% figure; 
% plot(Te_test,1e-6.*rate); 
% xlabel('T_e (eV)'); 
% ylabel('<\sigma v> (m^3/s)')


lnT = log(Te_eV);
lnn = log(ne_cm3/1e8);  %n tilde in fits

N = size(H4data,1);
M = size(H4data,2);

lnsv = zeros(size(Te_eV));
for n = 1:N
    for m = 1:M        
        lnsv = lnsv + H4data(m,n)*lnn^(n-1)*lnT.^(m-1);
    end
end

val = exp(lnsv);