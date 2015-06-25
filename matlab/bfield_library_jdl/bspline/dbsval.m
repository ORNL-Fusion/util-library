function sval = dbsval(x,kx,xknot,nx,bcoef)
% !
% !  Evaluates a spline, given its B-spline representation.
% !
% !   x      - point at which the spline is to be evaluated.  (input)
% !   kx     - order of the spline.  (input)
% !   xknot  - array of length nx+kx containing the knot
% !            sequence.  (input)
% !            xknot must be nondecreasing.
% !   nx     - number of B-spline coefficients.  (input)
% !   bcoef  - array of length nx containing the B-spline
% !            coefficients.  (input)
% !   dbsval - value of the spline at x.  (output)
% !

%     integer, intent(in)                          :: nx, kx
%     real(kind=dbl)                               :: dbsval
%     real(kind=dbl)                               :: x
%     real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
%     real(kind=dbl), dimension(nx), intent(in)    :: bcoef
% 
%     integer                       :: il, ik, ix, leftx
%     real(kind=dbl)                :: save1, save2
%     real(kind=dbl), dimension(kx) :: work, dl, dr

% !
% !     check if xknot(i) <= xknot(i+1) and calculation of i so that
% !     xknot(i) <= x < xknot(i+1)
% !

work = zeros(1,kx);
dl = zeros(1,kx);
dr = zeros(1,kx);

leftx = 0;

for ix = 1:nx+kx-1
    if (xknot(ix) > xknot(ix+1))
        fprintf('subroutine dbsval:\n')
        fprintf('xknot(ix) <= xknot(ix+1) required.\n')
        fprintf('ix = %i, xknot(ix) = %f, xknot(ix+1) = %f\n',ix,xknot(ix),xknot(ix+1))
        error('here')
    end
    if((xknot(ix) <= x) && (x < xknot(ix+1)))
        leftx = ix;
    end
end

if(leftx == 0)
    fprintf('subroutine dbsval:\n')
    fprintf('ix with xknot(ix) <= x < xknot(ix+1) required.')
    fprintf('x = %f\n', x)
    error('here')
end

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

sval = work(kx);

