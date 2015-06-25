function [work2,work3,bcoef] = spli2d(xyvec,ld,xydata,xyknot,n,k,m)
%     integer, intent(in)                         :: ld, n, k, m
%     real(kind=dbl), dimension(n), intent(in)    :: xyvec
%     real(kind=dbl), dimension(n+k), intent(in)  :: xyknot
%     real(kind=dbl), dimension(ld,m), intent(in) :: xydata
%     real(kind=dbl), dimension(m,n), intent(out) :: bcoef
%
%     real(kind=dbl), dimension(n), intent(out)         :: work2
%     real(kind=dbl), dimension((2*k-1)*n), intent(out) :: work3
%
%
%     integer        :: np1, km1, kpkm2, left, lenq, i, iflag, ilp1mx, j, jj
%     real(kind=dbl) :: xyveci

np1   = n + 1;
km1   = k - 1;
kpkm2 = 2 * km1;
left  = k;
lenq  = n * (k + km1);

work2 = zeros(1,n);
work3 = zeros(1,(2*k-1)*n);
bcoef = zeros(m,n);

for i = 1:lenq
    work3(i) = 0;
end

for i = 1:n
    xyveci  = xyvec(i);
    ilp1mx = min(i+k,np1);
    left   = max(left,i);
    if (xyveci < xyknot(left))
        fprintf('subroutine db2in:\n')
        fprintf('i with knot(i) <= x/y < knot(i+1) required.\n')
        fprintf('knot(1)   = %f\n', xyknot(1))
        fprintf('knot(n+k) = %f\n', xyknot(n+k))
        fprintf('      x/y = %f\n', xyveci)
        error('here')
    end
    while (xyveci >= xyknot(left+1)) 
        left = left + 1;
        if (left >= ilp1mx)
            left = left - 1;
            if (xyveci > xyknot(left+1))
                fprintf('subroutine db2in:\n')
                fprintf('i with knot(i) <= x/y < knot(i+1) required.\n')
                fprintf('knot(1)   = %f\n', xyknot(1))
                fprintf('knot(n+k) = %f\n', xyknot(n+k))
                fprintf('      x/y = %f\n', xyveci)
                error('here')
            end
        end
    end
    
    work2 = bsplvb(xyknot,n+k,k,1,xyveci,left);  %40
    jj = i - left + 1 + (left - k) * (k + km1);
    for j = 1:k
        jj        = jj + kpkm2;
        work3(jj) = work2(j);
    end
end

[work3,iflag] = banfac(work3,k+km1,n,km1,km1);

if (iflag ~= 1)
    fprintf('subroutine dbs2in: error\n')
    fprintf('no solution of linear equation system !!!\n')
    error('here')
end

for j = 1: m
    for i = 1: n
        work2(i) = xydata(i,j);
    end
    
    work2 = banslv(work3,k+km1,n,km1,km1,work2);
    
    for i = 1:n
        bcoef(j,i) = work2(i);
    end
end

