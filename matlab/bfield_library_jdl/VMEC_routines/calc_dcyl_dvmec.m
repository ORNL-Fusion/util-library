function [dRdu,dRdv,dZdu,dZdv] = calc_dcyl_dvmec(s,u,v,wout)

cosuv = cos(wout.xm*u-wout.xn*v);
sinuv = sin(wout.xm*u-wout.xn*v);

rmn_interp = ppval(wout.splrmn,s);
zmn_interp = ppval(wout.splzmn,s);

if wout.asym
    rmns_interp = ppval(wout.splrmns,s);
    zmnc_interp = ppval(wout.splzmnc,s);       
end

if nargout > 2
    npts = length(s);
    rs = rmn_interp.*sinuv;
    zc = zmn_interp.*cosuv;
    mr = repmat(wout.xm,1,npts);
    nr = repmat(wout.xn,1,npts);
    dRdu = sum(-mr.*rs);
    dRdv = sum(nr.*rs);
    dZdu = sum(mr.*zc);
    dZdv = sum(-nr.*zc);
    if wout.asym
        rc = rmns_interp.*cosuv;
        zs = zmnc_interp.*sinuv;        
        dRdu = dRdu + sum(mr.*rc);
        dRdv = dRdv + sum(-nr.*rc);
        dZdu = dZdu + sum(-mr.*zs);
        dZdv = dZdv + sum(nr.*zs);        
    end
end
