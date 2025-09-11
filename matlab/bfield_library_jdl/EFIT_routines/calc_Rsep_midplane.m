function [Rsep,Rsepi] = calc_Rsep_midplane(g)


Ninterp = 100;
Rinterp = linspace(g.rmaxis,g.r(end-2),Ninterp);
Zinterp = g.zmaxis*ones(1,Ninterp);
[psiN_mid,~] = calc_psiN(g,Rinterp,Zinterp);
Rsep = interp1(psiN_mid,Rinterp,1);  


Rinterp2 = linspace(g.rmaxis,g.r(2),Ninterp);
[psiN_midi,~] = calc_psiN(g,Rinterp2,Zinterp);
Rsepi = interp1(psiN_midi,Rinterp2,1);  

