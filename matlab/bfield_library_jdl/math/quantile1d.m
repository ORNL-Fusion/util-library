function q = quantile1d(x,p)
    n = length(x);    
    if nargin == 1
        p = linspace(0.5,n-0.5,n)./n;
    end
    xq = linspace(0.5,n-0.5,n)./n;     
    xs = sort(x);
    q = interp1(xq,xs,p);
    q(p<xq(1)) = xs(1);
    q(p>xq(end)) = xs(end);