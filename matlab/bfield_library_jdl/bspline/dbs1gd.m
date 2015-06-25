function val = dbs1gd(iderx,nxvec,xvec,kx,xknot,nx,bcoef)
%
% !
% !  Evaluates the derivative of a spline on a grid, given its B-spline
% !  representation.
% !
% !   iderx  - order of the derivative to be evaluated.  (input)
% !            in particular, iderx = 0 returns the value of the
% !            spline.
% !   nxvec  - length of vector xvec.  (input)
% !   xvec   - array of length nxvec containing the points at which the
% !            spline is to be evaluated.  (input)
% !            xvec should be strictly increasing.
% !   kx     - order of the spline.  (input)
% !   xknot  - array of length nx+kx containing the knot
% !            sequence.  (input)
% !            xknot must be nondecreasing.
% !   nx     - number of B-spline coefficients.  (input)
% !   bcoef  - array of length nx containing the B-spline
% !            coefficients.  (input)
% !   val    - array of length nxvec containing the values of the
% !            iderx-th derivative of the spline at the points in
% !            xvec.  (output)
% !

%     integer, intent(in)                           :: iderx, nxvec, kx, nx
%     real(kind=dbl), dimension(nxvec), intent(in)  :: xvec
%     real(kind=dbl), dimension(nx), intent(in)     :: bcoef
%     real(kind=dbl), dimension(nx+kx), intent(in)  :: xknot
%     real(kind=dbl), dimension(nxvec), intent(out) :: val
%
%     integer                             :: i, il, ik, ix
%     integer, dimension(nxvec)           :: leftx
%     real(kind=dbl)                      :: dik
%     real(kind=dbl), dimension(nxvec,kx) :: dl, dr, biatx, work
%     real(kind=dbl), dimension(nxvec)    :: save1, save2, term
%
%     logical :: same, next

leftx = zeros(1,nxvec);
biatx = zeros(nxvec,kx);
val = zeros(1,nxvec);
dr = zeros(nxvec,kx);
dl = zeros(nxvec,kx);
work = zeros(nxvec,kx);
save1 = zeros(1,nxvec);
save2 = zeros(1,nxvec);
term = zeros(1,nxvec);

leftx(1) =  huntn(xknot,nx+kx,kx,xvec(1),leftx(1));

for ix = 2:nxvec
    leftx(ix) = leftx(ix-1);
    same = (xknot(leftx(ix)) <= xvec(ix)) && (xvec(ix) <= xknot(leftx(ix)+1));
    if (~same )
        leftx(ix) = leftx(ix) + 1;
        next      = (xknot(leftx(ix)) <= xvec(ix)) && (xvec(ix) <= xknot(leftx(ix)+1));
        if (~next)
            leftx(ix) = huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix));
        end
    end
end

for ix = 1:nx+kx-1
    if (xknot(ix) > xknot(ix+1))
        fprintf('subroutine dbs1gd:\n')
        fprintf('xknot(ix) <= xknot(ix+1) required.\n')
        fprintf('ix=%f, xknot(ix)=%f, xknot(ix+1)=%f\n',ix, xknot(ix), xknot(ix+1))
        fprintf('\n')
        fprintf('xknot = %f\n',xknot)
        error('here')
    end
end

for ix = 1:nxvec
    if ((xvec(ix)<xknot(1))||(xvec(ix)>xknot(nx+kx)))
        fprintf('subroutine dbs1gd:\n')
        fprintf('ix with xknot(ix) <= x < xknot(ix+1) required.\n')
        fprintf('x = %f\n', xvec(ix))
        error('here')
    end
end

if (iderx == 0)
    
    for ix = 1:nxvec
        biatx(ix,1) = 1;
        val(ix)     = 0;
    end
    
    for ik = 1: kx-1
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
    
    for ik = 1:kx
        for ix = 1:nxvec
            val(ix) = val(ix) + biatx(ix,ik) * bcoef(leftx(ix)-kx+ik);
        end
    end
    
elseif ((iderx >= 1) && (iderx < kx))
    
    for ix = 1:nxvec
        biatx(ix,1) = 1;
        val(ix)     = 0;
    end
    
    for ik = 1:kx-iderx-1
        for ix = 1:nxvec
            dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix);
            dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+1-ik);
            save1(ix)    = biatx(ix,1);
            biatx(ix,1) = 0;
            for il = 1:ik
                term(ix)       = save1(ix)/ (dr(ix,il) + dl(ix,ik+1-il));
                biatx(ix,il)   = biatx(ix,il) + dr(ix,il) * term(ix);
                save1(ix)      = biatx(ix,il+1);
                biatx(ix,il+1) = dl(ix,ik+1-il) * term(ix);
            end
        end
    end
    
    for ik = 1:kx
        for ix = 1:nxvec
            work(ix,ik) = bcoef(leftx(ix)+ik-kx);
            dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix);
            dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+ik-kx);
        end
    end
    
    for ik = 1:iderx
        dik   = (kx - ik);
        for ix = 1:nxvec
            save2(ix) = work(ix,ik);
            for il = ik+1:kx
                save1(ix)   = work(ix,il);
                work(ix,il) = dik * (work(ix,il) - save2(ix))/(dl(ix,il) + dr(ix,il-ik));
                save2(ix)   = save1(ix);
            end
        end
    end
    
    for i = 1:kx-iderx
        for ix = 1:nxvec
            val(ix) = val(ix) + biatx(ix,i) * work(ix,iderx+i);
        end
    end
    
else
    
    for ix = 1:nxvec
        val(ix) = 0;
    end
    
end

