function val = dbs2gd(iderx,idery,nxvec,xvec,nyvec,yvec,kx,ky,xknot,yknot,nx,ny,bcoef,ldf)
% !
% !  Evaluates the derivative of a two-dimensional tensor-product spline,
% !  given its tensor-product B-spline representation on a grid.
% !
% !   iderx   - order of the derivative in the x-direction.  (input)
% !   idery   - order of the derivative in the y-direction.  (input)
% !   nxvec   - number of grid points in the x-direction.  (input)
% !   xvec    - array of length nx containing the x-coordinates at
% !             which the spline is to be evaluated.  (input)
% !             the points in xvec should be strictly increasing.
% !   nyvec   - number of grid points in the y-direction.  (input)
% !   yvec    - array of length ny containing the y-coordinates at
% !             which the spline is to be evaluated.  (input)
% !             the points in yvec should be strictly increasing.
% !   kx      - order of the spline in the x-direction.  (input)
% !   ky      - order of the spline in the y-direction.  (input)
% !   xknot   - array of length nx+kx containing the knot
% !             sequence in the x-direction.  (input)
% !             xknot must be nondecreasing.
% !   yknot   - array of length ny+ky containing the knot
% !             sequence in the y-direction.  (input)
% !             yknot must be nondecreasing.
% !   nx      - number of B-spline coefficients in the x-direction.
% !             (input)
% !   ny      - number of B-spline coefficients in the y-direction.
% !             (input)
% !   bcoef   - array of length nx*ny containing the
% !             tensor-product B-spline coefficients.  (input)
% !             bscoef is treated internally as a matrix of size nx
% !             by ny.
% !   val     - value of the (iderx,idery) derivative of the spline on
% !             the nx by ny grid.  (output)
% !             value(i,j) contains the derivative of the spline at the
% !             point (xvec(i),yvec(j)).
% !   ldf     - leading dimension of value exactly as specified in the
% !             dimension statement of the calling program.  (input)
% !

%     integer, intent(in)                           :: iderx, idery
%     integer, intent(in)                           :: nxvec, nyvec
%     integer, intent(in)                           :: kx, nx, ky, ny
%     integer, intent(in)                           :: ldf
%
%     real(kind=dbl), dimension(nxvec), intent(in)  :: xvec
%     real(kind=dbl), dimension(nyvec), intent(in)  :: yvec
%     real(kind=dbl), dimension(nx+kx), intent(in)  :: xknot
%     real(kind=dbl), dimension(ny+ky), intent(in)  :: yknot
%     real(kind=dbl), dimension(nx,ny), intent(in)  :: bcoef
%     real(kind=dbl), dimension(ldf,*), intent(out) :: val
%
%     integer                                     :: i, ik, il, ix, iy, ikx, iky
%     integer, dimension(nxvec)                   :: leftx
%     integer, dimension(nyvec)                   :: lefty
%     real(kind=dbl), dimension(nxvec,kx)         :: dl, dr
%     real(kind=dbl), dimension(max(nxvec,nyvec)) :: save1
%     real(kind=dbl), dimension(nxvec,kx)         :: biatx
%     real(kind=dbl), dimension(nyvec,ky)         :: biaty
%     real(kind=dbl), dimension(max(nxvec,nyvec)) :: term
%     real(kind=dbl), dimension(ky)               :: work
%
%     logical :: same,next


leftx = zeros(1,nxvec);
lefty = zeros(1,nyvec);
dl = zeros(nxvec,kx);
dr = zeros(nxvec,kx);
save1 = zeros(1,max(nxvec,nyvec));
biatx = zeros(nxvec,kx);
biaty = zeros(nyvec,ky);
term = zeros(1,max(nxvec,nyvec));
work = zeros(1,ky);
val = zeros(nxvec,nyvec);

leftx(1) = huntn(xknot,nx+kx,kx,xvec(1),leftx(1));

for ix = 2:nxvec
    leftx(ix) = leftx(ix-1);
    same = (xknot(leftx(ix)) <= xvec(ix)) && (xvec(ix) <= xknot(leftx(ix)+1));
    if(~same )
        leftx(ix) = leftx(ix) + 1;
        next      = (xknot(leftx(ix)) <= xvec(ix)) && (xvec(ix) <= xknot(leftx(ix)+1));
        if (~next)
            leftx(ix) = huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix));
        end
    end
end

for i = 1:nx+kx-1
    if (xknot(i) > xknot(i+1))
        fprintf('subroutine dbs2gd:\n')
        fprintf('xknot(i) <= xknot(i+1) required.\n')
        fprintf('i=%i, xknot(i)=%f, xknot(i+1)=%f\n',i, xknot(i), xknot(i+1))
        fprintf('\n')
        fprintf('xknot=%f\n',xknot)
        error('here')
    end
end

for i = 1:xvec
    if ((xvec(i)<xknot(1))||(xvec(i)>xknot(nx+kx)))
        fprintf('subroutine dbs2gd:\n')
        fprintf('ix with xknot(ix) <= x < xknot(ix+1) required.\n')
        fprintf('x = %f\n', xvec(i))
        error('here')
    end
end

lefty(1) = 0;

lefty(1)= huntn(yknot,ny+ky,ky,yvec(1),lefty(1));

for iy = 2:nyvec
    lefty(iy) = lefty(iy-1);
    same = (yknot(lefty(iy)) <= yvec(iy)) && (yvec(iy) <= yknot(lefty(iy)+1));
    if(~same )
        lefty(iy) = lefty(iy) + 1;
        next      = (yknot(lefty(iy)) <= yvec(iy))  && (yvec(iy) <= yknot(lefty(iy)+1));
        if (~next)
            lefty(iy)= huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy));
        end
    end
end

for i = 1:ny+ky-1
    if (yknot(i) > yknot(i+1))
        fprintf('subroutine dbs2gd:\n')
        fprintf('yknot(i) <= yknot(i+1) required.\n')
        fprintf('i=%i, yknot(i)=%f, yknot(i+1)=%f\n',i, yknot(i), yknot(i+1))
        fprintf('\n')
        fprintf('yknot=%f\n',yknot)
        error('here')
    end
end

for i = 1:nyvec
    if ((yvec(i)<yknot(1))||(yvec(i)>yknot(ny+ky)))
        fprintf('subroutine dbs2gd:\n')
        fprintf('iy with yknot(iy) <= y < yknot(iy+1) required.\n')
        fprintf('y = %f\n', yvec(i))
        error('here')
    end
end

if ((iderx == 0) && (idery == 0))
    
    for ix = 1:nxvec
        biatx(ix,1) = 1;
    end
    
    for ik = 1:kx-1
        for ix = 1:nxvec
            dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix);
            dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik);
            save1(ix) = 0;
        end
        
        for il = 1:ik
            for ix = 1:nxvec
                term(ix)     = biatx(ix,il)/ (dr(ix,il) + dl(ix,ik+1-il));
                biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix);
                save1(ix)    = dl(ix,ik+1-il) * term(ix);
            end
        end
        
        for ix = 1:nxvec
            biatx(ix,ik+1) = save1(ix);
        end
    end
    
    for iy = 1:nyvec
        biaty(iy,1) = 1;
    end
    
    for ik = 1:ky-1
        for iy = 1:nyvec
            dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy);
            dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik);
            save1(iy) = 0;
        end
        
        for il = 1:ik
            for iy = 1:nyvec
                term(iy)     = biaty(iy,il)/ (dr(iy,il) + dl(iy,ik+1-il));
                biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy);
                save1(iy)    = dl(iy,ik+1-il) * term(iy);
            end
        end
        
        for iy = 1:nyvec
            biaty(iy,ik+1) = save1(iy);
        end
    end
    
    for iy = 1:nyvec
        for ix = 1:nxvec
            val(ix,iy) = 0;
        end
    end
    
    for iky = 1:ky
        for ikx = 1:kx
            for iy = 1:nyvec
                for ix = 1:nxvec
                    val(ix,iy) = val(ix,iy)+ biatx(ix,ikx) * biaty(iy,iky)*bcoef(leftx(ix)-kx+ikx,lefty(iy)-ky+iky);
                end
            end
        end
    end
    
elseif (((iderx >= 1) || (idery >= 1)) && ( (iderx < kx) && (idery < ky)))
    
    for iy = 1:nyvec
        for ix = 1:nxvec
            for iky = 1:ky
                work(iky) = dbsdca(iderx,xvec(ix),kx,xknot,nx,bcoef(:,lefty(iy)-ky+iky),leftx(ix));
            end
            val(ix,iy) = dbsder(idery,yvec(iy),ky,yknot(lefty(iy)-ky+1:lefty(iy)+ky),ky,work);
        end
    end
    
else
    
    for iy = 1:nyvec
        for ix = 1:nxvec
            val(ix,iy) = 0;
        end
    end
    
end
