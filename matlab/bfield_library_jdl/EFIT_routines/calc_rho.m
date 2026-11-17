function [rho,Phi,Reff,Reff_bry] = calc_rho(g,R1,Z1)
%  [rho,Phi,Reff,Reff_bry] = calc_rho(g,R1,Z1)
% rho is sqrt normalized toroidal flux
% Phi is total toroidal flux (Wb)
% Reff satisfies pi*Reff^2 = abs(Phi/g.bcentr)

% See test_toroidal_flux 

%% Map to psiN first
psiN = calc_psiN(g,R1,Z1);
if any(psiN >= 1)
    if length(psiN) == 1
        error('rho undefined outside of separatrix. PsiN here = %f\n',psiN)
    else
        disp('Some values of rho requested outside of separatrix!')
    end
end

% Integrate q using poloidal flux per radian; 2*pi gives total toroidal flux.
PhiEval = g.ip_sign*cumtrapz(g.qpsi)*2*pi*(g.ssibry-g.ssimag)/(g.mw-1);
Phi = interp1(g.pn,PhiEval,psiN);
rho = sqrt(Phi/PhiEval(end));
Reff = sqrt(abs(Phi/(pi*g.bcentr)));
Reff_bry = sqrt(abs(PhiEval(end)/(pi*g.bcentr)));



