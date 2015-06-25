function dbs2dre
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !
% ! IMSL name:  dbs2dr (double precision version)
% !
% ! purpose:    evaluate the derivative of a two-dimensional
% !             tensor-product spline, given its tensor-product
% !             B-spline representation.
% !
% ! usage:      dbs2dr(ixder, iyder, x, y, kxord, kyord, xknot, yknot,
% !                    nxcoef, nycoef, bscoef)
% !
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !
% !        specifications for parameters
% !
kxord = 5;
kyord = 3;
nxdata = 21;
nydata = 6;
ldf = nxdata;

%       integer    i, j, nxcoef, nycoef
%
%       double precision bscoef(nxdata,nydata), f, f21,                   &
%      &     fdata(ldf,nydata), s21, x, xdata(nxdata),                    &
%      &     xknot(nxknot), y, ydata(nydata), yknot(nyknot)

% !
% !        define function and (2,1) derivative
% !

xdata = zeros(1,nxdata);
ydata = zeros(1,nydata);
fdata = zeros(ldf,nydata);
%
% !
% !        set up interpolation points
% !

for i = 1:nxdata
    xdata(i) = (i-11)/10.0;
end

% !
% !        generate knot sequence
% !

xknot = dbsnak (nxdata, xdata, kxord);

% !
% !        set up interpolation points
% !

for i = 1:nydata
    ydata(i) = (i-1)/5.0;
end

% !
% !        generate knot sequence
% !

yknot = dbsnak (nydata, ydata, kyord);

% !
% !        generate fdata
% !

for i = 1:nydata
    for j = 1:nxdata
        fdata(j,i) = f(xdata(j),ydata(i));
    end
end

% !
% !        interpolate
% !

bscoef = dbs2in (nxdata, xdata, nydata, ydata, fdata, ldf, kxord, kyord, xknot, yknot);

nxcoef = nxdata;
nycoef = nydata;

% !
% !        write heading
% !

fprintf('                              (2,1)\n')
fprintf('   x     y    (x,y)     error\n')

% !
% !        print (2,1) derivative over a grid of [0.0,1.0] x [0.0,1.0]
% !        at 16 points.
% !

for i = 1:4
    for j = 1:4
        x   = (i-1)/3.0;
        y   = (j-1)/3.0;
        
        % !
        % !              evaluate spline
        % !
        
        s21 = dbs2dr(2,1,x,y,kxord,kyord,xknot,yknot,nxcoef,nycoef,bscoef);
        fprintf('%15.4f   %15.4f  %15.4f  %15.6f\n', x, y, s21, f21(x,y) - s21)
    end
end


    function f = f(x,y)
        f= x*x*x*x + x*x*x*y*y;
    end
    function f21 = f21(x,y)
        f21 = 12.0*x*y;
    end

end
%
% !                                        (2,1)
% !              x              y          s    (x,y)     error
% !          0.0000         0.0000         0.0000       0.000000
% !          0.0000         0.3333         0.0000       0.000000
% !          0.0000         0.6667         0.0000       0.000000
% !          0.0000         1.0000         0.0000       0.000001
% !          0.3333         0.0000         0.0000       0.000000
% !          0.3333         0.3333         1.3333       0.000002
% !          0.3333         0.6667         2.6667      -0.000002
% !          0.3333         1.0000         4.0000       0.000008
% !          0.6667         0.0000         0.0000       0.000006
% !          0.6667         0.3333         2.6667      -0.000011
% !          0.6667         0.6667         5.3333       0.000028
% !          0.6667         1.0000         8.0001      -0.000134
% !          1.0000         0.0000        -0.0004       0.000439
% !          1.0000         0.3333         4.0003      -0.000319
% !          1.0000         0.6667         7.9996       0.000363
% !          1.0000         1.0000        12.0005      -0.000458


% MY FORTRAN OUTPUT::: JDL
%                                        (2,1)
%              x              y          s    (x,y)     error
%          0.0000         0.0000        -0.0000       0.000000
%          0.0000         0.3333         0.0000      -0.000000
%          0.0000         0.6667         0.0000      -0.000000
%          0.0000         1.0000         0.0000      -0.000000
%          0.3333         0.0000         0.0000       0.000000
%          0.3333         0.3333         1.3333      -0.000000
%          0.3333         0.6667         2.6667       0.000000
%          0.3333         1.0000         4.0000      -0.000000
%          0.6667         0.0000        -0.0000       0.000000
%          0.6667         0.3333         2.6667      -0.000000
%          0.6667         0.6667         5.3333      -0.000000
%          0.6667         1.0000         8.0000       0.000000
%          1.0000         0.0000         0.0000      -0.000000
%          1.0000         0.3333         4.0000      -0.000000
%          1.0000         0.6667         8.0000      -0.000000
%          1.0000         1.0000        12.0000      -0.000000
% 
