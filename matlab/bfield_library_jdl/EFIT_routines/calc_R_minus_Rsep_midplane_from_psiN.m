function [RminusRsep,R,Rsep] = calc_R_minus_Rsep_midplane_from_psiN(g,psiN)

nInterp = 100;
rInterp = linspace(g.rmaxis,g.r(end-1),nInterp);
zInterp = g.zmaxis*ones(size(rInterp));

psiNInterp = calc_psiN(g,rInterp,zInterp);

if psiN > max(psiNInterp) | psiN < min(psiNInterp)
    error('PsiN value out of range %f',psiN)
end
R = interp1(psiNInterp,rInterp,psiN);
Rsep = interp1(psiNInterp,rInterp,1);

RminusRsep = R - Rsep;

