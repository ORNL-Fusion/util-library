function f = gauss2d(x, y, A, x0, y0, sigx, sigy)
% GAUSS2D - Evaluates a 2D Gaussian function at given points (x, y)
%
% Inputs:
%   x    - x-coordinate(s) at which to evaluate the Gaussian (scalar or array)
%   y    - y-coordinate(s) at which to evaluate the Gaussian (scalar or array)
%   A    - Amplitude of the Gaussian (peak value)
%   x0   - x-coordinate of the center of the Gaussian
%   y0   - y-coordinate of the center of the Gaussian
%   sigx - Standard deviation (width) of the Gaussian in the x-direction
%   sigy - Standard deviation (width) of the Gaussian in the y-direction
%
% Output:
%   f    - Value of the Gaussian function at the input coordinates
%
% Formula:
%   f(x, y) = A * exp( -( (x - x0)^2 / (2 * sigx^2) + (y - y0)^2 / (2 * sigy^2) ) )
%
% Total area under curve:
%   V = 2 * pi * sigx * sigy
%   (Useful if normalization is needed)

    % Evaluate the 2D Gaussian
    f = A * exp( -( (x - x0).^2 ./ (2 * sigx^2) + (y - y0).^2 ./ (2 * sigy^2) ) );
end
