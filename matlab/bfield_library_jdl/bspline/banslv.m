function b = banslv(w,nroww,nrow,nbandl,nbandu,b)


%     integer, intent(in)                               :: nroww,nrow
%     integer, intent(in)                               :: nbandl,nbandu
%     real(kind=dbl), dimension(nroww,nrow), intent(in) :: w
%     real(kind=dbl), dimension(nrow), intent(inout)    :: b
% 
%     integer :: middle, nrowm1, jmax, i, j
w = reshape(w,nroww,nrow);

middle = nbandu + 1;
if (nrow == 1)
    b(1) = b(1) / w(middle,1);
    return;
end
nrowm1 = nrow - 1;
if (nbandl ~= 0)    
    for i = 1:nrowm1
        jmax = min(nbandl, nrow - i);
        for j = 1:jmax
            b(i+j) = b(i+j) - b(i) * w(middle+j,i);
        end
    end    
end
    
if (nbandu <= 0)
    for i = 1:nrow
       b(i) = b(i) / w(1,i);
    end
    return;
end

for i = nrow:-1:2
    b(i) = b(i)/w(middle,i);
    jmax = min(nbandu,i-1);
    for j = 1:jmax
        b(i-j) = b(i-j) - b(i) * w(middle-j,i);
    end
end
