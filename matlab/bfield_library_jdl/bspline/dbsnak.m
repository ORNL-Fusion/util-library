function xknot=dbsnak(nx,xvec,kxord)
%
%  Compute the `not-a-knot' spline knot sequence.
%  (see de Boor p. 167)
%
%   nx     - number of data points.  (input)
%   xvec   - array of length ndata containing the location of the
%            data points.  (input)
%   kxord  - order of the spline.  (input)
%   xknot  - array of length ndata+korder containing the knot
%            sequence.  (output)
%

% initialize
xknot = zeros(1,nx+kxord);
    
if kxord < 0 || kxord > nx
    fprintf('subroutine dbsnak: error\n')
    fprintf('0 <= kxord <= nx is required.\n')
    fprintf('kxord = %i, and nx = %i, is given.\n', kxord, nx)
    error('Error in dbsnak')
end

xknot(1:kxord) = xvec(1);

if mod(kxord,2) == 0
    for ix = kxord+1:nx
        xknot(ix) = xvec(ix-kxord/2);
    end
else
    for ix = kxord+1:nx
        xknot(ix) = 0.5*(xvec(ix-floor(kxord/2)) + xvec(ix-floor(kxord/2)-1));
    end
end

for ix = nx+1:nx+kxord
    xknot(ix) = xvec(nx) * (1 + eps);
end

