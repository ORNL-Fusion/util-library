function [Bx,By,Bz] = bfield_bs_jc( x,y,z,coil,current )

npts = size(coil,1);

nobs = length(x);
bx = zeros(nobs);
by = zeros(nobs);
bz = zeros(nobs);

c = current(1:npts-1)*1d-7;
for i = 1:nobs
    bl = -coil(1:npts-1,1)+x(i);
    bm = -coil(1:npts-1,2)+y(i);
    bn = -coil(1:npts-1,3)+z(i);
    bb = 1./sqrt(bl.^2 + bm.^2 + bn.^2);

    cl = -coil(2:npts,1)+x(i);
    cm = -coil(2:npts,2)+y(i);
    cn = -coil(2:npts,3)+z(i);
    cc = 1./sqrt(cl.^2 + cm.^2 + cn.^2);

    al = bl-cl;
    am = bm-cm;
    an = bn-cn;

    adotb = (al.*bl+am.*bm+an.*bn);
    adotc = (al.*cl+am.*cm+an.*cn);

    ul = cm.*bn-cn.*bm;
    um = cn.*bl-cl.*bn;
    un = cl.*bm-cm.*bl;

    recu = 1./(ul.^2 + um.^2 + un.^2);

    w = (adotc.*cc-adotb.*bb).*recu;   

    Bx(i) = sum(c.*w.*ul);
    By(i) = sum(c.*w.*um);
    Bz(i) = sum(c.*w.*un);   
end
end

