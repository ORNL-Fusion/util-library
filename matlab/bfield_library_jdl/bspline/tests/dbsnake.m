function dbsnake
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !
% ! IMSL name:  dbsnak (double precision version)
% !
% ! purpose:    compute the `not-a-knot' spline knot sequence.
% !
% ! usage:      call dbsnak(ndata, xdata, korder, xknot)
% !
% ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndata = 20;
kmin = 3;
kmax = 8;

DEBUG = 1;

xdata = zeros(1,ndata);
fdata = zeros(1,ndata);

for i = 1:ndata
    xt = (i-1)/(ndata-1);
    xdata(ndata-i+1) = t(xt);
end

xdata(1) = 0;

for i = 1:ndata
    fdata(i) = f(xdata(i));
end


if DEBUG
    figure; hold on;
    plot(xdata,fdata)
end
fprintf(' korder   max diff\n');
for k=kmin:kmax
    korder = k;
    % gen knots
    xknot = dbsnak(ndata,xdata,korder);
    % interpolate
    bscoef = dbsint(ndata, xdata, fdata, korder, xknot);
   
    difmax = 0;
    for i = 1:100
        xt = (i-1)/99.1;
        
        % evaluate spline
        st = dbsval(xt,korder,xknot,ndata,bscoef);
        ft = f(xt);
        dif = abs(ft-st);
        if DEBUG
            plot(xt,st,'ko')
        end
        % compute max diff
        
        difmax = max(dif,difmax);
        
    end
   
    fprintf('%i    %9.4f   \n',korder,difmax)
end
    
    function f=f(x)
        f=sin(10*x*x*x);
    end


    function t=t(x)
        t = 1-x*x;
    end
end


% !  korder     maximum difference
% !
% !    3        0.0081
% !    4        0.0026
% !    5        0.0004
% !    6        0.0008
% !    7        0.0010
% !    8        0.0004
