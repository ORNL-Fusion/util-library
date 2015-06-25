function bcoef = dbsint(nx,xvec,xdata,kx,xknot)
%
%  Computes the spline interpolant, returning the B-spline coefficients.
%  (see de Boor p. 204)
%
%   nx     - number of data points.  (input)
%   xvec   - array of length nx containing the data point
%            abscissas.  (input)
%   xdata  - array of length ndata containing the data point
%            ordinates.  (input)
%   kx     - order of the spline.  (input)
%            korder must be less than or equal to ndata.
%   xknot  - array of length nx+kx containing the knot
%            sequence.  (input)
%            xknot must be nondecreasing.
%   bscoef - array of length ndata containing the B-spline
%            coefficients.  (output)
%

bcoef = zeros(1,nx);
%
%     integer, intent(in)                          :: nx, kx
%     real(kind=dbl), dimension(nx), intent(in)    :: xdata, xvec
%     real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
%     real(kind=dbl), dimension(nx), intent(out)   :: bcoef
%
%     integer                                :: nxp1, kxm1, kpkm2, leftx, lenq
%     integer                                :: ix, ik,ilp1mx, jj, iflag
%     real(kind=dbl)                         :: xveci
%     real(kind=dbl), dimension((2*kx-1)*nx) :: work


nxp1  = nx + 1;
kxm1  = kx - 1;
kpkm2 = 2*kxm1;
leftx = kx;
lenq  = nx*(kx + kxm1);

work(1:lenq) = 0;

for  ix = 1:nx
    xveci  = xvec(ix);
    ilp1mx = min(ix+kx,nxp1);
    leftx   = max(leftx,ix);
    if (xveci < xknot(leftx))
        fprintf('subroutine dbsint:\n')
        fprintf('xknot(ix) <= xknot(ix+1) required.\n')
        fprintf('ix = %i, xknot(ix) = %f, xnot(ix+1) = %f\n',ix,xknot(ix),xknot(ix+1))
        error('here')
    end
    
    while (xveci >= xknot(leftx+1)) 
        leftx = leftx + 1;
        if (leftx >= ilp1mx)
            leftx = leftx - 1;
            if (xveci > xknot(leftx+1))
                fprintf('subroutine dbsint:\n')
                fprintf('xknot(ix) <= xknot(ix+1) required.\n')
                fprintf('ix = %i, xknot(ix) = %f, xnot(ix+1) = %f\n',ix,xknot(ix),xknot(ix+1))
                error('here')
            end
        end
    end
    
    bcoef= bsplvb(xknot,nx+kx,kx,1,xveci,leftx);
    jj = ix - leftx + 1 + (leftx - kx) * (kx + kxm1);
    for ik = 1:kx
        jj       = jj + kpkm2;
        work(jj) = bcoef(ik);
    end        
end

[work,iflag] = banfac(work,kx+kxm1,nx,kxm1,kxm1);

if (iflag ~= 1)
    fprintf('subroutine dbsint: error\n')
    fprintf('no solution of linear equation system !!!\n')
    error('Error in dbsint')
end

for ix = 1:nx
    bcoef(ix) = xdata(ix);
end

bcoef = banslv(work,kx+kxm1,nx,kxm1,kxm1,bcoef);