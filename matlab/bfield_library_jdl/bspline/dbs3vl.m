function sval = dbs3vl(x,y,z,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef)
% !
% !  Evaluates a three-dimensional tensor-product spline, given its
% !  tensor-product B-spline representation.
% !
% !   x      - x-coordinate of the point at which the spline is to be
% !            evaluated.  (input)
% !   y      - y-coordinate of the point at which the spline is to be
% !            evaluated.  (input)
% !   z      - z-coordinate of the point at which the spline is to be
% !            evaluated.  (input)
% !   kx     - order of the spline in the x-direction.  (input)
% !   ky     - order of the spline in the y-direction.  (input)
% !   kz     - order of the spline in the z-direction.  (input)
% !   xknot  - array of length nx+kx containing the knot
% !            sequence in the x-direction.  (input)
% !            xknot must be nondecreasing.
% !   yknot  - array of length ny+ky containing the knot
% !            sequence in the y-direction.  (input)
% !            yknot must be nondecreasing.
% !   zknot  - array of length nz+kz containing the knot
% !            sequence in the z-direction.  (input)
% !            zknot must be nondecreasing.
% !   nx     - number of B-spline coefficients in the x-direction.
% !            (input)
% !   ny     - number of B-spline coefficients in the y-direction.
% !            (input)
% !   nz     - number of B-spline coefficients in the z-direction.
% !            (input)
% !   bcoef  - array of length nx*ny*nz containing the
% !            tensor-product B-spline coefficients.  (input)
% !            bscoef is treated internally as a matrix of size nx
% !            by ny by nz.
% !   dbs3vl - value of the spline at (x,y,z).  (output)
% !

% !
% !     check if knot(i) <= knot(i+1) and calculation of i so that
% !     knot(i) <= x < knot(i+1)
% !

nintz = 0;
for iz = 1:nz+kz-1
    if (zknot(iz) > zknot(iz + 1))
        fprintf('subroutine dbs3vl:\n')
        fprintf('zknot(iz) <= zknot(iz+1) required.\n')
        fprintf('iz = %i, zknot(iz) = %f, zknot(iz+1) = %f\n',iz,zknot(iz),zknot(iz+1))
        error('here')
    end
    if ((zknot(iz) <= z) && (z < zknot(iz + 1)))
        nintz = iz;
    end
end

if (nintz == 0)
    fprintf('subroutine dbs3vl:\n')
    fprintf('iz with zknot(iz) <= z < zknot(iz+1) required.\n')
    fprintf('zknot(iz)   = %f\n', zknot(iz))
    fprintf('  z         = %f\n', z)
    fprintf('zknot(iz+1) = %f\n', zknot(iz+1))
    error('here')
end

for iz = 1:kz
    work(iz) = dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef(:,:,nintz-kz+iz));
end

sval = dbsval(z,kz,zknot(nintz-kz+1:nintz+kz),kz,work);

