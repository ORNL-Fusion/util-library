function [xshift,yshift] = find_IR_center(cells,data)
% Find the point that maximizes the summed data within a circle of radius
% radius_eval

use_fmincon = 1;
radius_eval = 9; % [cm]
dx_step = 0.25;  % [cm] - dx for gridded search of rel step size for fmincon
dx_minmax = 5;   % [cm]
debug = 0;
opt_display = 'off';     % Controls optimization output.  'off','final','iter'

if use_fmincon == 0
    % just a gridded search    
    xtest = -dx_minmax:dx_step:dx_minmax;
    f = zeros(length(xtest));
    for i = 1:length(xtest)
        for j = 1:length(xtest)
            f(i,j) = maxfun([xtest(i),xtest(j)]);
        end
    end    
    [mx,ix] = max(f);
    [~,iy] = max(mx);
    ix = ix(iy);
    xshift = xtest(ix);
    yshift = xtest(iy);
else    
    x00 = [0,0];
    lb = dx_minmax*[-1,-1];
    ub = dx_minmax*[1,1];
    opts=optimoptions(@fmincon,'findiffrelstep',dx_step,'Display',opt_display);
    xfinal = fmincon(@minfun,x00,[],[],[],[],lb,ub,[],opts);
    xshift = xfinal(1);
    yshift = xfinal(2);
end

    
    function f = minfun(shift)
        f = -maxfun(shift);
    end
    
    function f = maxfun(shift)
        xshift0 = shift(1);
        yshift0 = shift(2);
        rmesh = sqrt((cells.xmesh-xshift0).^2 + (cells.ymesh-yshift0).^2);
        f = sum(sum(data.IRdata2D(rmesh<radius_eval)));       
        if debug
            fprintf('x=%e, y=%e, f=%f\n',shift,f)
        end        
    end

end