function [r2,z2,pn2] = refine_psiN(r,z,fac,g)
%--------------------------------
%% Refine g psi
%--------------------------------
r2 = linspace(r(1),r(end),length(r)*fac);
z2 = linspace(z(1),z(end),length(z)*fac);
for i = 1:length(z2)
    pn2(:,i) = calc_psiN(g,r2,z2(i)*ones(size(r2)));
end
end