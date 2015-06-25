function dbsdere
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !
% ! IMSL name:  dbsder (double precision version)
% !
% ! purpose:    evaluate the derivative of a spline, given its B-spline
% !             representation.
% !
% ! usage:      dbsder(ideriv, x, korder, xknot, ncoef, bscoef)
% !
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !
% !        specifications for parameters
% !
korder = 3;
ndata = 5;

xdata = zeros(1,ndata);
fdata = zeros(1,ndata);

% !
% !       define function and derivative
% !
%
%       f(x)  = sqrt(x)
%       df(x) = 0.5/sqrt(x)

% !
% !        set up interpolation points
% !

for i = 1:ndata
    xdata(i) = i/ndata;
    fdata(i) = f(xdata(i));
end

% !
% !        generate knot sequence
% !

xknot = dbsnak (ndata, xdata, korder);

% !
% !        interpolate
% !

bscoef = dbsint (ndata, xdata, fdata, korder, xknot);

% !
% !        write heading
% !

fprintf('x    s(x)     error     s''(x)     error\n')

% !
% !        print on a finer grid
% !

ncoef = ndata;
xt    = xdata(1);

% !
% !        evaluate spline
% !

bt0   = dbsder(0,xt,korder,xknot,ncoef,bscoef);
bt1   = dbsder(1,xt,korder,xknot,ncoef,bscoef);

fprintf('  %6.4f   %7.4f    %10.6f    %8.4f     %10.6f\n',xt, bt0, f(xt) - bt0, bt1, df(xt) - bt1)

for i = 2:ndata
    xt  = (xdata(i-1)+xdata(i))/2.0;
    
    % !
    % !           evaluate spline
    % !
    
    bt0 = dbsder(0,xt,korder,xknot,ncoef,bscoef);
    bt1 = dbsder(1,xt,korder,xknot,ncoef,bscoef);
    fprintf('  %6.4f   %7.4f    %10.6f    %8.4f     %10.6f\n',xt, bt0, f(xt) - bt0, bt1, df(xt) - bt1)
    
    xt  = xdata(i);
    
    % !
    % !           evaluate spline
    % !
    
    bt0 = dbsder(0,xt,korder,xknot,ncoef,bscoef);
    bt1 = dbsder(1,xt,korder,xknot,ncoef,bscoef);
    fprintf('  %6.4f   %7.4f    %10.6f    %8.4f     %10.6f\n',xt, bt0, f(xt) - bt0, bt1, df(xt) - bt1)
end





    function f=f(x)
        f = sqrt(x);
    end
    function df=df(x)
        df = 0.5/sqrt(x);
    end
end
%
% !       x        s(x)       error        s'(x)        error
% !
% !  0.2000      0.4472     0.000000       1.0423     0.075738
% !  0.3000      0.5456     0.002084       0.9262    -0.013339
% !  0.4000      0.6325     0.000000       0.8101    -0.019553
% !  0.5000      0.7077    -0.000557       0.6940     0.013071
% !  0.6000      0.7746     0.000000       0.6446     0.000869
% !  0.7000      0.8366     0.000071       0.5952     0.002394
% !  0.8000      0.8944     0.000000       0.5615    -0.002525
% !  0.9000      0.9489    -0.000214       0.5279    -0.000818
% !  1.0000      1.0000     0.000000       0.4942     0.005814
