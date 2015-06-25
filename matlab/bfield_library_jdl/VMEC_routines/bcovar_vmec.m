function b = bcovar_vmec(s,u,v,wout)

% cosuv = cos(wout.xm*u-wout.xn*v);
cosuv = cos(wout.xm_nyq*u-wout.xn_nyq*v);

bsupu_interp = ppval(wout.splbsupu,s);
bsupv_interp = ppval(wout.splbsupv,s);

b.bsupu = sum(bsupu_interp.*cosuv);
b.bsupv = sum(bsupv_interp.*cosuv);

if wout.asym
    sinuv = sin(wout.xm_nyq*u-wout.xn_nyq*v);
    bsupus_interp = ppval(wout.splbsupus,s);
    bsupvs_interp = ppval(wout.splbsupvs,s);
    b.bsupu = b.bsupu + sum(bsupus_interp.*sinuv);
    b.bsupv = b.bsupv + sum(bsupvs_interp.*sinuv);
end
