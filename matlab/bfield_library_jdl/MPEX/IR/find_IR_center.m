function [xshift,yshift] = find_IR_center(cells,data)


    radius_eval = 9;
    dxtest = 0.5;
    xtest = -5:dxtest:5;
    for i = 1:length(xtest)
        for j = 1:length(xtest)             
            f(i,j) = maxfun(xtest(i),xtest(j));
        end
    end
    
    [mx,ix] = max(f);
    [my,iy] = max(mx);
    ix = ix(iy);
    xshift = xtest(ix);
    yshift = xtest(iy);
    function f = maxfun(xshift,yshift)
        rmesh = sqrt((cells.xmesh-xshift).^2 + (cells.ymesh-yshift).^2);
        f = sum(sum(data.IRdata2D(rmesh<radius_eval)));
    end

end