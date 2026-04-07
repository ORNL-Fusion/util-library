function [Br,Bz,Atheta] = bfield_circular_coils_vectorized(coil,current,r,z)

if ~isequal(size(r),size(z))
    error('r and z must have the same size')
end

nwind = length(coil.rwind);
npts = numel(r);

rvec = r(:);
zvec = z(:);
awind = reshape(coil.rwind,1,nwind);
dwind = reshape(coil.zwind,1,nwind);
cur = reshape(current,1,nwind);

rgrid = repmat(rvec,1,nwind);
zgrid = repmat(zvec,1,nwind);
agrid = repmat(awind,npts,1);
dgrid = repmat(dwind,npts,1);

mu0_4pi = 1.e-7;

Br_mat = zeros(npts,nwind);
Bz_mat = zeros(npts,nwind);
if nargout > 2
    Atheta_mat = zeros(npts,nwind);
end

is_axis = (rgrid == 0);
is_off_axis = ~is_axis;

if any(is_axis(:))
    dz0 = zgrid(is_axis) - dgrid(is_axis);
    a0 = agrid(is_axis);
    Bz_mat(is_axis) = 2*pi*mu0_4pi*a0.^2 ./ (a0.^2 + dz0.^2).^(3/2);
end

if any(is_off_axis(:))
    ro = rgrid(is_off_axis);
    zo = zgrid(is_off_axis);
    ao = agrid(is_off_axis);
    do = dgrid(is_off_axis);

    m = 4*ro.*ao ./ ((ro+ao).^2 + (zo-do).^2);
    sm = sqrt(m);
    [K,E] = ellipke(m);

    sa = sqrt(ao);
    sroa = sqrt(ro./ao);
    dz = zo - do;
    m1 = m - 1;

    Br_mat(is_off_axis) = -2*mu0_4pi./ro.^(1.5).*dz./sa.*(sm/2.*K - sm/4.*(m-2)./m1.*E);
    Bz_mat(is_off_axis) = 2*mu0_4pi./ro.*(sroa.*sm/2.*K - sm./(4*m1).*(sroa.*(m-2) + m./sroa).*E);

    if nargout > 2
        Atheta_mat(is_off_axis) = 2*mu0_4pi./sroa.*((2./sm - sm).*K - 2./sm.*E);
    end
end

Br = reshape(Br_mat*cur.',size(r));
Bz = reshape(Bz_mat*cur.',size(r));
if nargout > 2
    Atheta = reshape(Atheta_mat*cur.',size(r));
end
