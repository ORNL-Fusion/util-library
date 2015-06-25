function dbsinte
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !
% ! IMSL name:  dbsint (double precision version)
% !
% ! purpose:    compute the spline interpolant, returning the B-spline
% !             coefficients.
% !
% ! usage:      call dbsint(ndata, xdata, fdata, korder, xknot, bscoef)
% !
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !
% !        specifications for parameters
% !
korder = 3;
ndata = 5;

% !
% !        define function  (below)
% !
%       f(x) = sqrt(x)
%
% !
% !        set up interpolation points
% !

xdata = zeros(1,ndata);
fdata = zeros(1,ndata);

for i = 1:ndata
    xdata(i) = (i-1)/(ndata-1);
    fdata(i) = f(xdata(i));
end

% !
% !        generate knot sequence
% !

xknot = dbsnak(ndata, xdata, korder);

% !
% !       interpolate
% !

bscoef = dbsint (ndata, xdata, fdata, korder, xknot);

% !
% !        write heading
% !

fprintf('   x           s(x)           error\n');

% !
% !        print on a finer grid
% !

ncoef = ndata;
xt    = xdata(1);

% !
% !        evaluate spline
% !

bt = dbsval(xt,korder,xknot,ncoef,bscoef);

fprintf(' %6.4f       %8.4f    %11.6f\n', xt, bt, f(xt) - bt);

for i = 2:ndata
    xt = (xdata(i-1)+xdata(i))/2.0;
    
    % !
    % !           evaluate spline
    % !
    
    bt = dbsval(xt,korder,xknot,ncoef,bscoef);
    fprintf(' %6.4f       %8.4f    %11.6f\n', xt, bt, f(xt) - bt);
    
    
    xt = xdata(i);
    
    % !
    % !           evaluate spline
    % !
    
    bt = dbsval(xt,korder,xknot,ncoef,bscoef);
    fprintf(' %6.4f       %8.4f    %11.6f\n', xt, bt, f(xt) - bt);
    
end


    function f=f(x)
        f = sqrt(x);
    end

end
% !
% !       x                   s(x)                  error
% !
% !  0.0000                 0.0000               0.000000
% !  0.1250                 0.2918               0.061781
% !  0.2500                 0.5000               0.000000
% !  0.3750                 0.6247              -0.012311
% !  0.5000                 0.7071               0.000000
% !  0.6250                 0.7886               0.002013
% !  0.7500                 0.8660               0.000000
% !  0.8750                 0.9365              -0.001092
% !  1.0000                 1.0000               0.000000
