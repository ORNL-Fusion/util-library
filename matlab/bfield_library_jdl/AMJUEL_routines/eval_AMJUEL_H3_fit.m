function val = eval_AMJUEL_H3_fit(H4data,E0_eV,Ti_eV)
% Te in eV

lnT = log(Ti_eV);
lnE0 = log(E0_eV);  %n tilde in fits

N = size(H4data,1);
M = size(H4data,2);

lnsv = zeros(size(Ti_eV));
for n = 1:N
    for m = 1:M        
        lnsv = lnsv + H4data(m,n)*lnE0^(n-1)*lnT.^(m-1);
    end
end
% lnsv

val = exp(lnsv);