function dbs1gde
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !
% ! IMSL name:  dbs1gd (double precision version)
% !
% ! purpose:    evaluate the derivative of a spline on a grid, given
% !             its B-spline representation.
% !
% ! usage:      call dbs1gd(ideriv, n, xvec, korder, xknot, ncoef,
% !                         bscoef, value)
% !
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% !
% !        specifications for parameters
% !

korder = 3;
ndata = 5;
nknot = ndata + korder;
nfgrid = 9;

% !
% !        specifications for local variables
% !
%
%       integer    i, ncoef
%       double precision ans0(nfgrid), ans1(nfgrid), bscoef(ndata),       &
%      &           fdata(ndata),                                          &
%      &           x, xdata(ndata), xknot(nknot), xvec(nfgrid)
%


% !
% !        set up interpolation points
% !

xdata = zeros(1,ndata);
fdata = zeros(1,ndata);

for i = 1:ndata
    xdata(i) = (i)/(ndata);
    fdata(i) = f(xdata(i));
end

xknot =  dbsnak (ndata, xdata, korder);

% !
% !        interpolate
% !
bscoef =  dbsint (ndata, xdata, fdata, korder, xknot);
fprintf('   x         s(x)        error      s''(x)       error\n')

% !
% !        print on a finer grid
% !

ncoef   = ndata;
xvec(1) = xdata(1);

for i = 2:2:2*ndata - 2
    xvec(i)   = (xdata(i/2+1)+xdata(i/2))/2.0;
    xvec(i+1) = xdata(i/2+1);
end

ans0 = dbs1gd (0, 2*ndata-1, xvec, korder, xknot, ncoef, bscoef);
ans1 = dbs1gd (1, 2*ndata-1, xvec, korder, xknot, ncoef, bscoef);

for i = 1:2*ndata - 1
    fprintf('  %6.4f   %7.4f    %8.4f    %8.4f     %8.4f\n', xvec(i), ans0(i), f(xvec(i)) - ans0(i), ans1(i), df(xvec(i)) - ans1(i))
end


    function f=f(x)
        f = sqrt(x);
    end
    function df=df(x)
        df = 0.5/sqrt(x);
    end


end

% !       x        s(x)       error        s'(x)        error
% !
% !  0.2000      0.4472       0.0000       1.0423       0.0757
% !  0.3000      0.5456       0.0021       0.9262      -0.0133
% !  0.4000      0.6325       0.0000       0.8101      -0.0196
% !  0.5000      0.7077      -0.0006       0.6940       0.0131
% !  0.6000      0.7746       0.0000       0.6446       0.0009
% !  0.7000      0.8366       0.0001       0.5952       0.0024
% !  0.8000      0.8944       0.0000       0.5615      -0.0025
% !  0.9000      0.9489      -0.0002       0.5279      -0.0008
% !  1.0000      1.0000       0.0000       0.4942       0.0058
