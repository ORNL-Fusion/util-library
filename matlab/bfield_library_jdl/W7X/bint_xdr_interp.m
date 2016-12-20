function [bval,idiv]=bint_xdr_interp(bgrid,xvec)

r=xvec(1);
p=xvec(2);
z=xvec(3);

[poss,~,isign] = symmetrize([r,p,z],bgrid.nfp,bgrid.stellsym);
r = poss(1);
p = poss(2);
z = poss(3);

if r < bgrid.rmin || r > bgrid.rmax || z < bgrid.zmin || z > bgrid.zmax
    idiv = 1;
    bval(1:3) = NaN;
else
    idiv = 0;
    method = 'linear';
    error('This is too inefficient, find indices since grid is regular')
    br = interp3(bgrid.rg,bgrid.pg,bgrid.zg,bgrid.brg,r,p,z,method,0.);
    bp = interp3(bgrid.rg,bgrid.pg,bgrid.zg,bgrid.bfg,r,p,z,method,0.);
    bz = interp3(bgrid.rg,bgrid.pg,bgrid.zg,bgrid.bzg,r,p,z,method,0.);
%     bval = [vecrot*isign*[br,bp,bz].'].';
    bval = [isign*br,bp,bz];
%     bval = [[br,bp,bz].'].';
end
