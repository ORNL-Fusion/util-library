function [work2,work3,bcoef] = spli3d(xyzvec,ldf,mdf,xyzdata,xyzknot,n,k,m,l,nx,ny,nz)
%
%     integer, intent(in)                               :: ldf, mdf, n, k, m, l
%     integer, intent(in)                               :: nx, ny, nz
%     real(kind=dbl), dimension(n), intent(in)          :: xyzvec
%     real(kind=dbl), dimension(n+k), intent(in)        :: xyzknot
%     real(kind=dbl), dimension(ldf,mdf,*), intent(in)  :: xyzdata
%     real(kind=dbl), dimension(nx,ny,nz), intent(out)  :: bcoef
%     real(kind=dbl), dimension(n), intent(out)         :: work2
%     real(kind=dbl), dimension((2*k-1)*n), intent(out) :: work3
%
%     integer        :: np1, km1, kpkm2, left, lenq, i, ilp1mx, j, jj, iflag, in
%     real(kind=dbl) :: xyzveci

bcoef = zeros(nx,ny,nz);
work3 = zeros(1,(2*k-1)*n);

np1   = n + 1;
km1   = k - 1;
kpkm2 = 2 * km1;
left  = k;
lenq  = n * (k + km1);

work3(1:lenq) = 0;

for i = 1:n
    xyzveci = xyzvec(i);
    ilp1mx  = min(i+k,np1);
    left    = max(left,i);
    if (xyzveci < xyzknot(left))
        fprintf('subroutine spli3d:\n')
        fprintf('i with knot(i) <= x/y/z < knot(i+1) required.\n')
        fprintf('knot(1)   = %f\n', xyzknot(1))
        fprintf('knot(n+k) = %f\n', xyzknot(n+k))
        fprintf('    x/y/z = %f\n', xyzveci)
        error('here')
    end
    while (xyzveci >= xyzknot(left+1))
        left = left + 1;
        if (left >= ilp1mx)
            left = left - 1;
            if (xyzveci > xyzknot(left+1))
                fprintf('subroutine spli3d:\n')
                fprintf('i with knot(i) <= x/y/z < knot(i+1) required.\n')
                fprintf('knot(1)   = %f\n', xyzknot(1))
                fprintf('knot(n+k) = %f\n', xyzknot(n+k))
                fprintf('    x/y/z = %f\n', xyzveci)
                error('here')
            end
        end
    end
    work2 = bsplvb(xyzknot,n+k,k,1,xyzveci,left);
    jj = i - left + 1 + (left - k) * (k + km1);
    for j = 1:k
        jj    = jj + kpkm2;
        work3(jj) = work2(j);
    end
end

[work3,iflag] =  banfac(work3,k+km1,n,km1,km1);

if (iflag ~= 1)
    fprintf('subroutine dbs3in: error\n')
    fprintf('no solution of linear equation system !!!\n')
    error('here')
end

for j = 1:l
    for i = 1:m
        for in = 1:n
            work2(in) = xyzdata(i,j,in);
        end
        
        work2 = banslv(work3,k+km1,n,km1,km1,work2);
        
        for in = 1:n
            bcoef(i,j,in) = work2(in);
        end
        
    end
end




