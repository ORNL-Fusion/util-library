function dbs3dr = dbs3dr(iderx,idery,iderz,x,y,z,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef)
%
% !
% !  Evaluates the derivative of a three-dimensional tensor-product spline,
% !  given its tensor-product B-spline representation.
% !
% !   iderx  - order of the x-derivative.  (input)
% !   idery  - order of the y-derivative.  (input)
% !   iderz  - order of the z-derivative.  (input)
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
% !   dbs3dr - value of the (iderx,idery,iderz) derivative of the
% !            spline at (x,y,z).  (output)
% !

%     integer, intent(in)                              :: iderx, idery, iderz
%     integer, intent(in)                              :: nx, ny, nz, kx, ky, kz
%     real(kind=dbl), intent(in)                       :: x, y, z
%     real(kind=dbl), dimension(nx+kx), intent(in)     :: xknot
%     real(kind=dbl), dimension(ny+ky), intent(in)     :: yknot
%     real(kind=dbl), dimension(nz+kz), intent(in)     :: zknot
%     real(kind=dbl), dimension(nx,ny,nz), intent(in)  :: bcoef
%     real(kind=dbl)                                   :: dbs3dr
%
%     integer                       :: iz, nintz
%     real(kind=dbl), dimension(kz) :: work

% !
% !     check if knot(i) <= knot(i+1) and calculation of i so that
% !     knot(i) <= x < knot(i+1)
% !
work = zeros(1,kz);

nintz = 0;

for iz = 1:nz+kz-1
    if (zknot(iz) > zknot(iz + 1))
        fprintf('subroutine dbs3vr:\n')
        fprintf('zknot(iz) <= zknot(iz+1) required.\n')
        fprintf('iz=%i, zknot(iz)=%f, zknot(iz+1)=%f\n', iz, zknot(iz), zknot(iz+1))
        error('here')
    end
    if((zknot(iz) <= z) && (z < zknot(iz + 1)))
        nintz = iz;
    end
end

if(nintz == 0)
    fprintf('subroutine dbs3dr:\n')
    fprintf('iz with zknot(iz) <= z < zknot(iz+1) required.\n')
    fprintf('zknot(iz)   = %f\n', zknot(iz))
    fprintf('  z         = %f\n', z)
    fprintf('zknot(iz+1) = %f\n', zknot(iz+1))
    error('here')
end

for iz = 1:kz
    work(iz) = dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef(:,:,nintz-kz+iz));
end

dbs3dr = dbsder(iderz,z,kz,zknot(nintz-kz+1:nintz+kz),kz,work);
% 
% !                                       (2,0,1)
% !          x           y           z    s     (x,y,z)    error
% !      -1.0000      0.0000      0.0000   -0.000107    0.000107
% !      -1.0000      0.0000      1.0000    0.000053   -0.000053
% !      -1.0000      0.3333      0.0000    0.064051   -0.064051
% !      -1.0000      0.3333      1.0000   -5.935941   -0.064059
% !      -1.0000      0.6667      0.0000    0.127542   -0.127542
% !      -1.0000      0.6667      1.0000  -11.873034   -0.126966
% !      -1.0000      1.0000      0.0000    0.191166   -0.191166
% !      -1.0000      1.0000      1.0000  -17.808527   -0.191473
% !      -0.3333      0.0000      0.0000   -0.000002    0.000002
% !      -0.3333      0.0000      1.0000    0.000000    0.000000
% !      -0.3333      0.3333      0.0000    0.021228   -0.021228
% !      -0.3333      0.3333      1.0000   -1.978768   -0.021232
% !      -0.3333      0.6667      0.0000    0.042464   -0.042464
% !      -0.3333      0.6667      1.0000   -3.957536   -0.042464
% !      -0.3333      1.0000      0.0000    0.063700   -0.063700
% !      -0.3333      1.0000      1.0000   -5.936305   -0.063694
% !       0.3333      0.0000      0.0000   -0.000003    0.000003
% !       0.3333      0.0000      1.0000    0.000000    0.000000
% !       0.3333      0.3333      0.0000   -0.021229    0.021229
% !       0.3333      0.3333      1.0000    1.978763    0.021238
% !       0.3333      0.6667      0.0000   -0.042465    0.042465
% !       0.3333      0.6667      1.0000    3.957539    0.042462
% !       0.3333      1.0000      0.0000   -0.063700    0.063700
% !       0.3333      1.0000      1.0000    5.936304    0.063697
% !       1.0000      0.0000      0.0000   -0.000098    0.000098
% !       1.0000      0.0000      1.0000    0.000053   -0.000053
% !       1.0000      0.3333      0.0000   -0.063855    0.063855
% !       1.0000      0.3333      1.0000    5.936146    0.063854
% !       1.0000      0.6667      0.0000   -0.127631    0.127631
% !       1.0000      0.6667      1.0000   11.873067    0.126933
% !       1.0000      1.0000      0.0000   -0.191442    0.191442
% !       1.0000      1.0000      1.0000   17.807940    0.192060
