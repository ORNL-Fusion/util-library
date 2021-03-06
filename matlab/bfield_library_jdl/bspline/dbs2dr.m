function dbs2dr = dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef)
%
% !
% !  Evaluates the derivative of a two-dimensional tensor-product spline,
% !  given its tensor-product B-spline representation.
% !
% !   iderx  - order of the derivative in the x-direction.  (input)
% !   idery  - order of the derivative in the y-direction.  (input)
% !   x      - x-coordinate of the point at which the spline is to be
% !            evaluated.  (input)
% !   y      - y-coordinate of the point at which the spline is to be
% !            evaluated.  (input)
% !   kx     - order of the spline in the x-direction.  (input)
% !   ky     - order of the spline in the y-direction.  (input)
% !   xknot  - array of length nx+kx containing the knot
% !            sequence in the x-direction.  (input)
% !            xknot must be nondecreasing.
% !   yknot  - array of length ny+ky containing the knot
% !            sequence in the y-direction.  (input)
% !            yknot must be nondecreasing.
% !   nx     - number of B-spline coefficients in the x-direction.
% !            (input)
% !   ny     - number of B-spline coefficients in the y-direction.
% !            (input)
% !   bcoef  - array of length nx*ny containing the
% !            tensor-product B-spline coefficients.  (input)
% !            bscoef is treated internally as a matrix of size nx
% !            by ny.
% !   dbs2dr  - value of the (iderx,idery) derivative of the spline at
% !            (x,y).  (output)
% !
%
%     use numeric
%
%     implicit none
%
%     integer, intent(in)                          :: iderx, idery
%     integer, intent(in)                          :: kx, nx, ky, ny
%     real(kind=dbl)                               :: dbs2dr
%     real(kind=dbl), intent(in)                   :: x, y
%     real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
%     real(kind=dbl), dimension(ny+ky), intent(in) :: yknot
%     real(kind=dbl), dimension(nx,ny), intent(in) :: bcoef
%
%     integer                       :: ix, iy, iky, nintx, ninty
%     real(kind=dbl), dimension(ky) :: work
%
% !
% !     check if knot(i) <= knot(i+1) and calculation of i so that
% !     knot(i) <= x < knot(i+1)
% !

work = zeros(1,ky);

% nintx = 0;
% for ix = 1:nx+kx-1
%     if (xknot(ix) > xknot(ix+1))
%         fprintf('subroutine dbs2dr:\n')
%         fprintf('xknot(ix) <= xknot(ix+1) required.\n')
%         fprintf('ix=%i, xknot(ix)=%f, xknot(ix+1)=%f\n', ix, xknot(ix), xknot(ix+1))
%         error('here')
%     end
%     if((xknot(ix) <= x) && (x < xknot(ix+1)))
%         nintx = ix;
%     end
% end

if ~all(diff(xknot)>=0)
    test = diff(xknot);    
    ix = test(test<0);
    ix = ix(1);
    fprintf('subroutine dbs2dr:\n')
    fprintf('xknot(ix) <= xknot(ix+1) required.\n')
    fprintf('ix=%i, xknot(ix)=%f, xknot(ix+1)=%f\n', ix, xknot(ix), xknot(ix+1))
    error('here')
end
[nintx] = find((xknot(1:nx+kx-1) <= x) & (x < xknot(2:nx+kx)));

if(nintx == 0)
    fprintf('subroutine dbs2dr:\n')
    fprintf('ix with xknot(ix) <= x < xknot(ix+1) required.\n')
    fprintf('x = %f\n', x)
    error('here')
end

% ninty = 0;
% for iy = 1:ny+ky-1
%     if (yknot(iy) > yknot(iy+1))
%         fprintf('subroutine dbs2dr:\n')
%         fprintf('yknot(iy) <= yknot(iy+1) required.\n')
%         fprintf('iy=%i, yknot(iy)=%f, yknot(iy+1)=%f\n', iy, yknot(iy), yknot(iy+1))
%         error('here')
%     end
%     if ((yknot(iy) <= y) && (y < yknot(iy+1)))
%         ninty = iy;
%     end
% end

if ~all(diff(yknot)>=0)
    test = diff(yknot);    
    iy = test(test<0);
    iy = iy(1);
    fprintf('subroutine dbs2dr:\n')
    fprintf('yknot(iy) <= yknot(iy+1) required.\n')
    fprintf('iy=%i, yknot(iy)=%f, yknot(iy+1)=%f\n', iy, yknot(iy), yknot(iy+1))
    error('here')
end
[ninty] = find((yknot(1:ny+ky-1) <= y) & (y < yknot(2:ny+ky)));

if(ninty == 0)
    fprintf('subroutine dbs2dr:\n')
    fprintf('iy with yknot(iy) <= y < yknot(iy+1) required.\n')
    fprintf('y = %f\n', y)
    error('here')
end

for iky = 1: ky
    work(iky) =  dbsdca(iderx,x,kx,xknot,nx,bcoef(1:nx,ninty-ky+iky),nintx);
end

dbs2dr = dbsder(idery,y,ky,yknot(ninty-ky+1:ninty+ky),ky,work);
