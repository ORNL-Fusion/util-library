function f = eval_interp_test(x,y,testnum)
% Defines Franke [1979] 2D test functions (to be evaluated x,y = [0,1])
switch testnum
    case 1
        % Combination of four Gaussians
        f =   0.75*exp(-((9*x - 2).^2     + (9*y - 2).^2)/4 ) ...
            + 0.75*exp(-((9*x + 1).^2)/49 - (9*y + 1)/10)     ...
            + 0.50*exp(-((9*x - 7).^2     + (9*y - 3).^2)/4)  ...
            - 0.20*exp(-(9*x - 4).^2 - (9*y - 7).^2);
    case 2
        % hyperbolic tangent (cliff)
        f = (tanh(9*y - 9*x) + 1)/9;
    case 3
        % Saddle shape
        f = (1.25 + cos(5.4*y))./(6*(1 + (3*x - 1).^2));
    case 4
        % Gaussian with gentle slope
        f = exp(-81*((x - 0.5).^2 + (y - 0.5).^2)/16)/3;
    case 5
        % Gaussian with steep slope
        f = exp(-81*((x - 0.5).^2 + (y - 0.5).^2)/4)/3;
    case 6
        % Part of a sphere
        f = sqrt(64 - 81*((x - 0.5).^2 + (y - 0.5).^2))/9 - 0.5;
    otherwise
        error('Bad value for testnum: %i',testnum)
end

