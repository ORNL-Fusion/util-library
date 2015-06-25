function [R,Z,dRdu,dRdv,dZdu,dZdv] = suv_to_rz(s,u,v,wout)

cosuv = cos(wout.xm*u-wout.xn*v);
sinuv = sin(wout.xm*u-wout.xn*v);

rmn_interp = ppval(wout.splrmn,s);
zmn_interp = ppval(wout.splzmn,s);

R = sum(rmn_interp.*cosuv);
Z = sum(zmn_interp.*sinuv);

if wout.asym
    rmns_interp = ppval(wout.splrmns,s);
    zmnc_interp = ppval(wout.splzmnc,s);    
    R = R + sum(rmns_interp.*sinuv);
    Z = Z + sum(zmnc_interp.*cosuv);    
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