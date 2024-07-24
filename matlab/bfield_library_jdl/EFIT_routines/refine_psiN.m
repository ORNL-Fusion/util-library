function [r2,z2,pn2] = refine_psiN(r,z,fac,g)
% [r2,z2,pn2] = refine_psiN(r,z,fac,g)

if fac <= 0
    error('Bad value for fac in refine_psiN, must be >= 0',fac)
end
%--------------------------------
%% Refine g psi
%--------------------------------
nr2 = length(r)*fac;
nz2 = length(z)*fac;
r2 = linspace(r(1),r(end),nr2);
z2 = linspace(z(1),z(end),nz2);

pn2 = nan(nr2,nz2);
for i = 1:length(z2)
    pn2(:,i) = calc_psiN(g,r2,z2(i)*ones(size(r2)));
end

end