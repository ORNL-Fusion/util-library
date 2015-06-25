function dbs2gde
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !
% ! IMSL name:  dbs2gd (double precision version)
% !
% ! purpose:    evaluate the derivative of a two-dimensional
% !             tensor-product spline, given its tensor-product
% !             B-spline representation on a grid.
% !
% ! usage:      call dbs2gd(ixder, iyder, nx, xvec, ny, yvec, kxord,
% !                         kyord, xknot, yknot, nxcoef, nycoef,
% !                         bscoef, value, ldvalu)
% !
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% !
% !        specifications for local variables
% !
%
%       integer    i, j, kxord, kyord, ldf, nxcoef, nxdata,               &
%      &           nycoef, nydata
%
%       double precision dccfd(21,6), docbsc(21,6), docxd(21), docxk(26), &
%      &           docyd(6), docyk(9), f, f21, value(4,4),                &
%      &           x, xvec(4), y, yvec(4)
%
% !
% !        define function and derivative
% !

%
% !
% !       initialize/setup
% !

kxord  = 5;
kyord  = 3;
nxdata = 21;
nydata = 6;
ldf    = nxdata;

% !
% !       set up interpolation points
% !

for i = 1: nxdata
    docxd(i) = (i-11)/10.0;
end

% !
% !       set up interpolation points
% !

for i = 1:nydata
    docyd(i) = (i-1)/5.0;
end

% !
% !       generate knot sequence
% !

docxk = dbsnak (nxdata, docxd, kxord);

% !
% !       generate knot sequence
% !

docyk =  dbsnak (nydata, docyd, kyord);

% !
% !       generate fdata
% !

for i = 1:nydata
    for j = 1:nxdata
        dccfd(j,i) = f(docxd(j),docyd(i));
    end
end

% !
% !        interpolate
% !

docbsc = dbs2in (nxdata, docxd, nydata, docyd, dccfd, ldf, kxord,kyord, docxk, docyk);

% !
% !       print (2,1) derivative over a grid of [0.0,1.0] x [0.0,1.0]
% !       at 16 points.
% !

nxcoef = nxdata;
nycoef = nydata;
fprintf('                  (2,1)\n');
fprintf('   x       y     s     (x,y)      error  \n');

for i = 1:4
    xvec(i) = (i-1)/3.0;
    yvec(i) = xvec(i);
end

value = dbs2gd (2, 1,4, xvec, 4, yvec, kxord, kyord, docxk, docyk, nxcoef, nycoef, docbsc, 4);
% val = dbs2gd (iderx,idery,nxvec,xvec,nyvec,yvec,kx,ky,xknot,yknot,nx,ny,bcoef,ldf)
for i = 1:4
    for j = 1:4
        fprintf('%13.4f %13.4f %13.4f %15.6f\n', xvec(i), yvec(j), value(i,j), f21(xvec(i),yvec(j)) - value(i,j));
    end
end

    function f = f(x,y)
        f   = x*x*x*x + x*x*x*y*y;
    end
    function f21 = f21(x,y)
        f21 = 12.0*x*y;
    end
end

% !                                        (2,1)
% !              x              y          s    (x,y)     error
% !          0.0000         0.0000         0.0000       0.000000
% !          0.0000         0.3333         0.0000       0.000000
% !          0.0000         0.6667         0.0000       0.000000
% !          0.0000         1.0000         0.0000       0.000001
% !          0.3333         0.0000         0.0000      -0.000001
% !          0.3333         0.3333         1.3333       0.000001
% !          0.3333         0.6667         2.6667      -0.000003
% !          0.3333         1.0000         4.0000       0.000008
% !          0.6667         0.0000         0.0000      -0.000001
% !          0.6667         0.3333         2.6667      -0.000009
% !          0.6667         0.6667         5.3333       0.000037
% !          0.6667         1.0000         8.0001      -0.000120
% !          1.0000         0.0000        -0.0005       0.000488
% !          1.0000         0.3333         4.0003      -0.000320
% !          1.0000         0.6667         7.9994       0.000610
% !          1.0000         1.0000        12.0005      -0.000488
