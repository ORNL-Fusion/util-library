function dbsder = dbsder(iderx,x,kx,xknot,nx,bcoef)
% !
% !  Evaluates the derivative of a spline, given its B-spline representation.
% !
% !
% !   iderx  - order of the derivative to be evaluated.  (input)
% !            in particular, iderx = 0 returns the value of the
% !            spline.
% !   x      - point at which the spline is to be evaluated.  (input)
% !   kx     - order of the spline.  (input)
% !   xknot  - array of length nx+kx containing the knot
% !            sequence.  (input)
% !            xknot must be nondecreasing.
% !   nx     - number of B-spline coefficients.  (input)
% !   bcoef  - array of length nx containing the B-spline
% !            coefficients.  (input)
% !   dbsder - value of the iderx-th derivative of the spline at x.
% !            (output)
% !

%     integer, intent(in)                          :: iderx, kx, nx
%     real(kind=dbl)                               :: dbsder
%     real(kind=dbl), intent(in)                   :: x
%     real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
%     real(kind=dbl), dimension(nx), intent(in)    :: bcoef
%
%     integer                       :: ix, ik, il, leftx
%     real(kind=dbl)                :: save, save1, save2, y, sum0, dik
%     real(kind=dbl), dimension(kx) :: work, dl, dr,bsp
%
% !
% !     check if xknot(i) <= xknot(i+1) and calculation of i so that
% !     xknot(i) <= x < xknot(i+1)
% !

work = zeros(1,kx);
dl = zeros(1,kx);
dr = zeros(1,kx);
bsp = zeros(1,kx);

leftx = 0;
for ix = 1:nx+kx-1
    if (xknot(ix) > xknot(ix+1))
        fprintf('subroutine dbsder:\n')
        fprintf('xknot(ix) <= xknot(ix+1) required.\n')
        error('here')
    end
    if ((xknot(ix) <= x) && (x < xknot(ix+1)))
        leftx = ix;
    end
end

if (leftx == 0)
    fprintf('subroutine dbsder:\n')
    fprintf('ix with xknot(ix) <= x < xknot(ix+1) required.\n')
    fprintf('xknot(1)     = %f\n', xknot(1))
    fprintf('xknot(nx+kx) = %f\n', xknot(nx+kx))
    fprintf('         x   = %f\n', x)
    error('here')
end

if (iderx == 0)
    for ik = 1:kx-1
        work(ik) = bcoef(leftx+ik-kx);
        dl(ik)   = x - xknot(leftx+ik-kx);
        dr(ik)   = xknot(leftx+ik) - x;
    end
    
    work(kx)  = bcoef(leftx);
    dl(kx)    = x - xknot(leftx);
    
    for ik = 1:kx-1
        save2 = work(ik);
        for il = ik+1:kx
            save1 = work(il);
            work(il) = (dl(il) * work(il) + dr(il-ik) * save2)/ (dl(il) + dr(il - ik));
            save2 = save1;
        end
    end
    
    dbsder = work(kx);
    
elseif ((iderx >= 1) && (iderx <= kx))
    
    bsp(1) = 1;
    for ik = 1:kx-iderx-1
        dr(ik) = xknot(leftx+ik) - x;
        dl(ik) = x - xknot(leftx+1-ik);
        save   = bsp(1);
        bsp(1) = 0;
        for il = 1:ik
            y         = save / (dr(il) + dl(ik+1-il));
            bsp(il)   = bsp(il) + dr(il) * y;
            save      = bsp(il+1);
            bsp(il+1) = dl(ik+1-il) * y;
        end
    end
    
    for ik = 1:kx
        work(ik) = bcoef(leftx+ik-kx);
        dr(ik)   = xknot(leftx+ik) - x;
        dl(ik)   = x - xknot(leftx+ik-kx);
    end
    
    for ik = 1:iderx
        dik   = (kx - ik);
        save2 = work(ik);
        for il = ik+1:kx
            save1    = work(il);
            work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik));
            save2    = save1;
        end
    end
    
    sum0 = 0;
    
    for ix = 1:kx-iderx
        sum0 = sum0 + bsp(ix) * work(iderx+ix);
    end
    
    dbsder = sum0;
    
else
    dbsder = 0;
end

