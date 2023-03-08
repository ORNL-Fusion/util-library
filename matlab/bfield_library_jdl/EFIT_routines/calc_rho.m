function [rho,Phi,Reff,Reff_bry] = calc_rho(g,R1,Z1)
%  [rho,Phi,Reff,Reff_bry] = calc_rho(g,R1,Z1)
% Rho is sqrt toroidal flux
% Phi is toroidal flux
% Reff is R from Bo*pi*R^2 = 2*pi*Phi

% See test_toroidal_flux 

%% Map to psiN first
psiN = calc_psiN(g,R1,Z1);
if psiN >= 1
    error('rho undefined outside of separatrix. PsiN here = %f\n',psiN)
end

% Integrate q = Phi'/Psi' and interpolate at psiN
PhiEval = g.ip_sign*cumtrapz(g.qpsi)*2*pi*(g.ssibry-g.ssimag)/(g.mw-1);
Phi = interp1(g.pn,PhiEval,psiN);
rho = sqrt(Phi/PhiEval(end));
Reff = sqrt(abs(2*pi*Phi/(pi*g.bcentr)));
Reff_bry = sqrt(abs(2*pi*PhiEval(end)/(pi*g.bcentr)));



