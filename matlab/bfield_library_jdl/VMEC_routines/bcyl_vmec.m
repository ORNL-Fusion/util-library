function [Br,Bz,Bphi] = bcyl_vmec(R,Z,phi,wout)

[s,u,v,dRds,dZds,dRdu,dRdv,dZdu,dZdv]=rzp_to_suv(R,Z,phi,wout);
b = bcovar_vmec(s,u,v,wout);

Br = dRdu*b.bsupu + dRdv*b.bsupv;
Bz = dZdu*b.bsupu + dZdv*b.bsupv;
Bphi = R*b.bsupv;

