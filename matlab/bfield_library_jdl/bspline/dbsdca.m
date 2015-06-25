function dbsdca = dbsdca(iderx,x,kx,xknot,nx,bcoef,leftx)
% !
% ! This routine is equivalent to the routine dbsder, but it does not
% ! check the parameters!!!
% !
% ! Evaluates the derivative of a spline, given its B-spline representation.
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
% !   leftx  - number of the intervall of xknot that includes x
% !   dbsdca - value of the ideriv-th derivative of the spline at x.
% !            (output)
% !


work = zeros(1,kx);
dl = zeros(1,kx);
dr = zeros(1,kx);
bsp = zeros(1,kx);

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
    
    dbsdca = work(kx);
    
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
    
    dbsdca = sum0;
    
else
    dbsdca = 0;
end

