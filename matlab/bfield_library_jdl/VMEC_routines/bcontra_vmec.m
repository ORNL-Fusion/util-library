function b = bcontra_vmec(s,u,v,wout)

cosuv = cos(wout.xm_nyq*u-wout.xn_nyq*v);
sinuv = sin(wout.xm_nyq*u-wout.xn_nyq*v);

bsubu_interp = ppval(wout.splbsubu,s);
bsubv_interp = ppval(wout.splbsubv,s);
bsubs_interp = ppval(wout.splbsubs,s);

b.bsubu = sum(bsubu_interp.*cosuv);
b.bsubv = sum(bsubv_interp.*cosuv);
b.bsubs = sum(bsubs_interp.*sinuv);

if wout.asym
    bsubus_interp = ppval(wout.splbsubus,s);
    bsubvs_interp = ppval(wout.splbsubvs,s);
    bsubc_interp = ppval(wout.splbsubc,s);    
    b.bsubu = b.bsubu + sum(bsubu_interp.*sinuv);
    b.bsubv = b.bsubv + sum(bsubv_interp.*sinuv);
    b.bsubs = b.bsubs + sum(bsubs_interp.*cosuv);
end

%     mr = repmat(wout.xm,1,npts);
%     nr = repmat(wout.xn,1,npts);
