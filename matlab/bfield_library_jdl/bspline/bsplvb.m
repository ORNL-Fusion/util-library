function biatx= bsplvb(t,n,jhigh,index,x,left)
%     integer, intent(in) :: n, jhigh, index, left
%
%     real(kind=dbl), intent(in)                    :: x
%     real(kind=dbl), dimension(n), intent(in)      :: t
%     real(kind=dbl), dimension(jhigh), intent(out) :: biatx
%
%     integer                          :: j = 1
%     integer                          :: i, jp1
%     real(kind=dbl)                   :: saved, term
%     real(kind=dbl), dimension(jhigh) :: dl, dr


biatx = zeros(1,jhigh);
dl = zeros(1,jhigh);
dr = zeros(1,jhigh);

if (index == 1)
    j = 1;
    biatx(1) = 1;
    if (j >= jhigh)
        return
    end
end

while j < jhigh
    jp1 = j + 1;
    
    dr(j) = t(left+j) - x;
    dl(j) = x - t(left+1-j);
    saved = 0;
    
    for i = 1:j
        term     = biatx(i) / (dr(i) + dl(jp1-i));
        biatx(i) = saved + dr(i) * term;
        saved    = dl(jp1-i) * term;
    end
    
    biatx(jp1) = saved;
    j          = jp1;
end

