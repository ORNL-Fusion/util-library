function [w,iflag] = banfac(w,nroww,nrow,nbandl,nbandu)
%     integer, intent(in)                                  :: nroww,nrow
%     integer, intent(in)                                  :: nbandl,nbandu
%     integer, intent(out)                                 :: iflag
%     real(kind=dbl), dimension(nroww,nrow), intent(inout) :: w
%
%     real(kind=dbl) :: pivot, factor
%     integer        :: middle, nrowm1, jmax, kmax, ipk, midmk, i, j, k

w = reshape(w,nroww,nrow);

iflag  = 1;
middle = nbandu + 1;
nrowm1 = nrow - 1;

if (nrowm1 < 0)
    iflag = 2;
    w = reshape(w,1,nroww*nrow);
    return;
end

if (nrowm1 == 0)
     if (w(middle,nrow) == 0)   
         iflag = 2;
     end
     w = reshape(w,1,nroww*nrow);
     return;
end

if nbandl <= 0
    for i = 1:nrowm1
        if w(middle,i) == 0
            iflag = 2;
            w = reshape(w,1,nroww*nrow);
            return;
        end
    end
    
    if (nrowm1 == 0)
        if (w(middle,nrow) == 0)
            iflag = 2;
        end
        w = reshape(w,1,nroww*nrow);
        return;
    end
end
if (nbandu <= 0)    
    for i = 1:nrowm1
        pivot = w(middle,i);
        if(pivot == 0)
            iflag = 2;
            w = reshape(w,1,nroww*nrow);
            return;
        end
        jmax = min(nbandl, nrow - i);
        for j = 1:jmax
            w(middle+j,i) = w(middle+j,i) / pivot;
        end
    end    
    return;
end

for i = 1:nrowm1
    pivot = w(middle,i);
    if (pivot == 0)
        iflag = 2;
        w = reshape(w,1,nroww*nrow);
        return;
    end
    jmax = min(nbandl,nrow - i);
    for j = 1:jmax
        w(middle+j,i) = w(middle+j,i) / pivot;
    end
                
    kmax = min(nbandu,nrow - i);
                
    for k = 1:kmax
        ipk    = i + k;
        midmk  = middle - k;
        factor = w(midmk,ipk);
        for j = 1:jmax
            w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor;
        end
    end
end
w = reshape(w,1,nroww*nrow);