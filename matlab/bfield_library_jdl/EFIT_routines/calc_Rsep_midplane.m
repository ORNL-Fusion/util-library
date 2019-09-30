function [Rsep] = calc_Rsep_midplane(g)


Ninterp = 100;
Rinterp = linspace(g.rmaxis,g.r(end-1),Ninterp);
Zinterp = g.zmaxis*ones(1,Ninterp);
[psiN_mid,~] = calc_psiN(g,Rinterp,Zinterp);

Rsep = interp1(psiN_mid,Rinterp,1);   
