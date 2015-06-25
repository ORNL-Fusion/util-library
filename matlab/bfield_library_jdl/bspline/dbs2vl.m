function sval = dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef)
                
% !
% !  evaluates a two-dimensional tensor-product spline, given its
% !  tensor-product B-spline representation.    use numeric
% !
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
% !   dbs2vl - value of the spline at (x,y).  (output)
% !

% !
% !     check if knot(i) <= knot(i+1) and calculation of i so that
% !     knot(i) <= x < knot(i+1)
% !
leftx = 0;
for ix = 1:nx+kx-1
    if (xknot(ix) > xknot(ix+1))
        fprintf('subroutine dbs2vl:\n')
        fprintf('xknot(ix) <= xknot(ix+1) required.\n')
        fprintf('ix = %i, xknot(ix) = %f, xknot(ix+1) = %f\n',ix,xknot(ix),xknot(ix+1))
        error('here')
    end
    if((xknot(ix) <= x) && (x < xknot(ix+1))) 
        leftx = ix;
    end
end

if(leftx == 0)
    fprintf('subroutine dbs2vl:\n')
    fprintf('ix with xknot(ix) <= x < xknot(ix+1) required.\n')
    fprintf('x = %f\n',x)
    error('here')
end
       
lefty = 0;

for iy = 1:ny+ky-1
    if (yknot(iy) > yknot(iy+1))
        fprintf('subroutine dbs2vl:\n')
        fprintf('yknot(iy) <= yknot(iy+1) required.\n')
        fprintf('iy = %i, yknot(iy) = %f, yknot(iy+1) = %f\n',iy,yknot(iy),yknot(iy+1))
        error('here')
    end
    if ((yknot(iy) <= y) && (y < yknot(iy+1)))
        lefty = iy;
    end
end
    
if (lefty == 0)
    fprintf('subroutine dbs2vl:\n')
    fprintf('iy with yknot(iy) <= y < yknot(iy+1) required.\n')
    fprintf('y = %f\n',y)
    error('here')
end
    
    
for iky = 1:ky
    work(iky) = dbsdca(0,x,kx,xknot,nx,bcoef(1:nx,lefty-ky+iky),leftx);
end

sval = dbsval(y,ky,yknot(lefty-ky+1:lefty+ky),ky,work);
