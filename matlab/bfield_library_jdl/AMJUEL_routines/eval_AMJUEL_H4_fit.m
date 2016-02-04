function val = eval_AMJUEL_H4_fit(H4data,ne_cm3,Te_eV)
% Te in eV

lnT = log(Te_eV);
lnn = log(ne_cm3/1e8);  %n tilde in fits

N = size(H4data,2);
M = size(H4data,2);

lnsv = zeros(size(Te_eV));
for n = 1:N
    for m = 1:M        
        lnsv = lnsv + H4data(m,n)*lnn^(n-1)*lnT.^(m-1);
    end
end
% lnsv

val = exp(lnsv);