function [Br,Bz] = bfield_circular_coils_vectorized_zconst(coil,current,r,z)
% Optimized for the common MPEX fieldline-following path where z is scalar
% and only Br,Bz are needed.

if ~isscalar(z)
    error('z must be scalar')
end

nwind = length(coil.rwind);
rvec = r(:);
awind = reshape(coil.rwind,1,nwind);
dwind = reshape(coil.zwind,1,nwind);
cur = reshape(current,1,nwind);

mu0_4pi = 1.e-7;

Br_sum = zeros(size(rvec));
Bz_sum = zeros(size(rvec));

is_axis = (rvec == 0);
is_off_axis = ~is_axis;

if any(is_axis(:))
    dz0 = z - dwind;
    Bz_axis = 2*pi*mu0_4pi*awind.^2 ./ (awind.^2 + dz0.^2).^(3/2);
    Bz_sum(is_axis) = Bz_axis*cur.';
end

if any(is_off_axis(:))
    ro = rvec(is_off_axis);
    dz = z - dwind;

    m = 4*(ro.*awind) ./ ((ro+awind).^2 + dz.^2);
    sm = sqrt(m);
    [K,E] = ellipke(m);

    sa = sqrt(awind);
    sroa = sqrt(ro./awind);
    m1 = m - 1;

    Br_terms = -2*mu0_4pi./ro.^(1.5).*dz./sa.*(sm/2.*K - sm/4.*(m-2)./m1.*E);
    Bz_terms = 2*mu0_4pi./ro.*(sroa.*sm/2.*K - sm./(4*m1).*(sroa.*(m-2) + m./sroa).*E);
    Br_sum(is_off_axis) = Br_terms*cur.';
    Bz_sum(is_off_axis) = Bz_terms*cur.';
end

Br = reshape(Br_sum,size(r));
Bz = reshape(Bz_sum,size(r));
