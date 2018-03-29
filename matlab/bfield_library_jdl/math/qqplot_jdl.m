function qqplot_jdl(x,y)
    nx = length(x);
    ny = length(y);
    n = min([nx,ny]);
    
    p = linspace(0.5,n-0.5,n)./n;
    
    qx = quantile1d(x,p);
    qy = quantile1d(y,p);
    
    figure; hold on; box on;
    plot(qx,qy,'x')
    xlabel('x quantiles')
    ylabel('y quantiles')
    
    