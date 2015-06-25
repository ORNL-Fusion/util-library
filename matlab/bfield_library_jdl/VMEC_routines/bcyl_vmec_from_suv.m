function [Br,Bz,Bphi] = bcyl_vmec_from_suv(s,u,v,wout)

b = bcovar_vmec(s,u,v,wout);
[R,Z,dRdu,dRdv,dZdu,dZdv] = suv_to_rz(s,u,v,wout);

Br = dRdu*b.bsupu + dRdv*b.bsupv;
Bz = dZdu*b.bsupu + dZdv*b.bsupv;
Bphi = R*b.bsupv;

