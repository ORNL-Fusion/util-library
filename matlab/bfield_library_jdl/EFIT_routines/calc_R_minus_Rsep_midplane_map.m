function [RminusRsep,R,Rsep] = calc_R_minus_Rsep_midplane_map(g,R1,Z1)
% Give R,Z 
% RminusRsep = R - Rsep

Ninterp = 100;
Rinterp = linspace(g.rmaxis,g.r(end-1),Ninterp);
Zinterp = g.zmaxis*ones(1,Ninterp);

[psiN_mid,psi_mid] = calc_psiN(g,Rinterp,Zinterp);
[psiN1,psi1] = calc_psiN(g,R1,Z1);

R_midplane_map = interp1(psiN_mid,Rinterp,psiN1);   

Rsep = interp1(psiN_mid,Rinterp,1);   
rho_midplane_map = R_midplane_map - Rsep;

RminusRsep = rho_midplane_map;
R = R_midplane_map;
