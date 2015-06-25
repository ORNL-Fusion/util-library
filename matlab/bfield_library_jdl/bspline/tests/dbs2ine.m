function dbs2ine
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !
% ! IMSL name:  dbs2in (double precision version)
% !
% ! purpose:    compute a two-dimensional tensor-product spline
% !             interpolant, returning the tensor-product B-spline
% !             coefficients.
% !
% ! usage:      call dbs2in(nxdata, xdata, nydata, ydata, fdata, ldf,
% !                         kxord, kyord, xknot, yknot, bscoef)
% !
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !
% !        specifications for parameters
% !
kxord = 5;
kyord = 2;
nxdata = 21;
nxvec = 4;
nydata = 6;
nyvec = 4;
ldf = nxdata;
% nxknot = nxdata + kxord;
% nyknot = nydata + kyord;

%       integer    i, j, nxcoef, nycoef
%
%       double precision bscoef(nxdata,nydata), f, fdata(ldf,nydata),     &
%      &           value(nxvec,nyvec), x, xdata(nxdata), xknot(nxknot),   &
%      &           xvec(nxvec), y, ydata(nydata), yknot(nyknot),          &
%      &           yvec(nyvec)
%
% !
% !        define function
% !

xvec = zeros(1,nxvec);
yvec = zeros(1,nyvec);
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

xknot =  dbsnak (nxdata, xdata, kxord);
%
% !
% !        set up interpolation points
% !

for i = 1:nydata
    ydata(i) = (i-1)/5.0;
end
%
% !
% !        generate knot sequence
% !

yknot = dbsnak (nydata, ydata, kyord);

% !
% !        generate fdata
% !

for i = 1: nydata
    for j = 1:nxdata
        fdata(j,i) = f(xdata(j),ydata(i));
    end
end

% !
% !        interpolate
% !

bscoef =  dbs2in (nxdata, xdata, nydata, ydata, fdata, ldf, kxord, kyord, xknot, yknot);

nxcoef = nxdata;
nycoef = nydata;

% !
% !                                  write heading
% !

fprintf('   x     y      s(x,y)    error\n')


% !
% !        print over a grid of
% !        [0.0,1.0] x [0.0,1.0] at 16 points.
% !

for i = 1: nxvec
    xvec(i) = (i-1)/3.0;
end

for i = 1: nyvec
    yvec(i) = (i-1)/3.0;
end

% !
% !        evaluate spline
% !

value = dbs2gd (0, 0, nxvec, xvec, nyvec, yvec, kxord, kyord, xknot,yknot, nxcoef, nycoef, bscoef, nxvec);

for i = 1:nxvec
    for j = 1:nyvec
        fprintf('%15.4f   %15.4f   %15.4f  %15.6f\n',xvec(i), yvec(j),value(i,j),f(xvec(i),yvec(j))-value(i,j))
    end
end



    function f=f(x,y)
        f = x*x*x + x*y;
    end

end
%
%
% !              x              y          s(x,y)         error
% !          0.0000         0.0000         0.0000       0.000000
% !          0.0000         0.3333         0.0000       0.000000
% !          0.0000         0.6667         0.0000       0.000000
% !          0.0000         1.0000         0.0000       0.000000
% !          0.3333         0.0000         0.0370       0.000000
% !          0.3333         0.3333         0.1481       0.000000
% !          0.3333         0.6667         0.2593       0.000000
% !          0.3333         1.0000         0.3704       0.000000
% !          0.6667         0.0000         0.2963       0.000000
% !          0.6667         0.3333         0.5185       0.000000
% !          0.6667         0.6667         0.7407       0.000000
% !          0.6667         1.0000         0.9630       0.000000
% !          1.0000         0.0000         1.0000       0.000000
% !          1.0000         0.3333         1.3333       0.000000
% !          1.0000         0.6667         1.6667       0.000000
% !          1.0000         1.0000         2.0000       0.000000
