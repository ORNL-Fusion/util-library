function dbs3gde
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% IMSL name:  dbs3gd (double precision version)
%
% purpose:    evaluate the derivative of a three-dimensional
%             tensor-product spline, given its tensor-product
%             B-spline representation on a grid.
%
% usage:      call dbs3gd(ixder, iyder, izder, nx, xvec, ny, yvec, nz,
%                         zvec, kxord, kyord, kzord, xknot, yknot,
%                         zknot, nxcoef, nycoef, nzcoef, bscoef, value,
%                         ldvalu, mdvalu)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       use bspline
%
%       implicit none
%
% %
% %        specifications for parameters
% %
%
%       integer    kxord, kyord, kzord, ldf, ldval, mdf, mdval, nxdata,   &
%      &           nxknot, nydata, nyknot, nz, nzdata, nzknot
%       parameter  (kxord=5, kyord=2, kzord=3, ldval=4, mdval=4,          &
%      &           nxdata=21, nydata=6, nz=2, nzdata=8, ldf=nxdata,       &
%      &           mdf=nydata, nxknot=nxdata+kxord, nyknot=nydata+kyord,  &
%      &           nzknot=nzdata+kzord)
%
%       integer    i, j, k, l, nxcoef, nycoef, nzcoef
%       double precision bscoef(nxdata,nydata,nzdata), f, f201,           &
%      &           fdata(ldf,mdf,nzdata), value(ldval,mdval,nz),          &
%      &           x, xdata(nxdata), xknot(nxknot), xvec(ldval), y,       &
%      &           ydata(nydata), yknot(nyknot), yvec(mdval), z,          &
%      &           zdata(nzdata), zknot(nzknot), zvec(nz)
%
% %
% %        define function.
% %

kxord=5; kyord=2; kzord=3; ldval=4; mdval=4;
nxdata=21; nydata=6; nz=2; nzdata=8; ldf=nxdata;
mdf=nydata; nxknot=nxdata+kxord; nyknot=nydata+kyord;
nzknot=nzdata+kzord;


%
%        set up x interpolation points
%
for i = 1:nxdata
    xdata(i) = 2.0*((i-1)/(nxdata-1)) - 1.0;
end

%
%        set up y interpolation points
%

for i = 1:nydata
    ydata(i) = (i-1)/(nydata-1);
end

%
%        set up z interpolation points
%

for i = 1:nzdata
    zdata(i) = (i-1)/(nzdata-1);
end

%
%        generate knots
%

xknot = dbsnak (nxdata, xdata, kxord);
yknot = dbsnak (nydata, ydata, kyord);
zknot = dbsnak (nzdata, zdata, kzord);

%
%        generate fdata
%

for k = 1:nzdata
    for i = 1:nydata
        for j = 1:nxdata
            fdata(j,i,k) = f(xdata(j),ydata(i),zdata(k));
        end
    end
end

%
%        interpolate
%

bscoef = dbs3in (nxdata, xdata, nydata, ydata, nzdata, zdata, fdata, ldf, mdf, kxord, kyord, kzord, xknot, yknot, zknot);

nxcoef = nxdata;
nycoef = nydata;
nzcoef = nzdata;

%
%        print over a grid of [-1.0,1.0] x [0.0,1.0] x [0.0,1.0]
%        at 32 points.
%

for i = 1:4
    xvec(i) = 2.0*((i-1)/3.0) - 1.0;
end

for j = 1:4
    yvec(j) = (j-1)/3.0;
end

for l = 1:2
    zvec(l) = (l-1);
end

value = dbs3gd (2, 0, 1, 4, xvec, 4, yvec, 2, zvec, kxord, kyord, kzord, xknot, yknot, zknot, nxcoef, nycoef, nzcoef, bscoef, ldval, mdval);

fprintf('                (2,0,1)\n')
fprintf('    x        y          z        s      (x,y,z)       error\n')

for i = 1:4
    for j = 1:4
        for l = 1:2
            fprintf('%13.4f %13.4f %13.4f %13.4f %13.4f\n', xvec(i), yvec(j), zvec(l), value(i,j,l), f201(xvec(i),yvec(j),zvec(l)) - value(i,j,l));
        end
    end
end


    function f = f(x,y,z)
        f    = x*x*x*x + x*x*x*y*z*z*z;
    end
    function f201 = f201(x,y,z)
        f201 = 18.0*x*y*z;
    end

end

%                                             (2,0,1)
%           x           y          z          s     (x,y,z)  error
%
%      -1.0000       0.0000       0.0000       0.0000       0.0000
%      -1.0000       0.0000       1.0000       0.0000       0.0000
%      -1.0000       0.3333       0.0000       0.0637      -0.0637
%      -1.0000       0.3333       1.0000      -5.9363      -0.0637
%      -1.0000       0.6667       0.0000       0.1274      -0.1274
%      -1.0000       0.6667       1.0000     -11.8726      -0.1274
%      -1.0000       1.0000       0.0000       0.1911      -0.1911
%      -1.0000       1.0000       1.0000     -17.8089      -0.1911
%      -0.3333       0.0000       0.0000       0.0000       0.0000
%      -0.3333       0.0000       1.0000       0.0000       0.0000
%      -0.3333       0.3333       0.0000       0.0212      -0.0212
%      -0.3333       0.3333       1.0000      -1.9788      -0.0212
%      -0.3333       0.6667       0.0000       0.0425      -0.0425
%      -0.3333       0.6667       1.0000      -3.9575      -0.0425
%      -0.3333       1.0000       0.0000       0.0637      -0.0637
%      -0.3333       1.0000       1.0000      -5.9363      -0.0637
%       0.3333       0.0000       0.0000       0.0000       0.0000
%       0.3333       0.0000       1.0000       0.0000       0.0000
%       0.3333       0.3333       0.0000      -0.0212       0.0212
%       0.3333       0.3333       1.0000       1.9788       0.0212
%       0.3333       0.6667       0.0000      -0.0425       0.0425
%       0.3333       0.6667       1.0000       3.9575       0.0425
%       0.3333       1.0000       0.0000      -0.0637       0.0637
%       0.3333       1.0000       1.0000       5.9363       0.0637
%       1.0000       0.0000       0.0000       0.0000       0.0000
%       1.0000       0.0000       1.0000       0.0000       0.0000
%       1.0000       0.3333       0.0000      -0.0637       0.0637
%       1.0000       0.3333       1.0000       5.9363       0.0637
%       1.0000       0.6667       0.0000      -0.1274       0.1274
%       1.0000       0.6667       1.0000      11.8726       0.1274
%       1.0000       1.0000       0.0000      -0.1911       0.1911
%       1.0000       1.0000       1.0000      17.8089       0.1911
