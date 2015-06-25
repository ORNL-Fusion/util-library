function jlo = huntn(xx,n,kord,x,jlo)

%     integer, intent(in)                      :: n, kord
%     real(kind=dbl), intent(in)               :: x
%     real(kind=dbl), dimension(n), intent(in) :: xx
%
%     integer, intent(out)                     :: jlo
%
%     integer :: max, null, jhi, jm, inc
%
% !
% !     works only for B-Splines (order n)
% !

max  = n - kord;
null = kord;

if (jlo <= null || jlo > max)
    jlo = null;
    jhi = max+1;
    while jhi-jlo ~= 1
        jm = floor((jhi + jlo) / 2);
        if (x > xx(jm))
            jlo = jm;
        else
            jhi = jm;
        end
    end
    return;
end

inc = 1;

if (x >= xx(jlo))
    jhi = jlo + inc;
    if (jhi > max)
        jhi = max + 1;
    elseif (x >= xx(jhi))
        while (x < xx(jhi))
            jlo = jhi;
            inc = inc + inc;
            jhi = jlo + inc;
            if (jhi > max)
                break;
            end
        end
    end
else
    jhi = jlo;
    jlo = jhi - inc;
    if (jlo <= null)
        jlo = null;
    elseif (x < xx(jlo))
        while x >= xx(jlo)
            jhi = jlo;
            inc = inc + inc;
            if (jlo <= null)
                jlo = null;
                break;
            end
        end
    end
end

while jhi-jlo ~= 1
    jm = floor((jhi + jlo) / 2);
    if (x > xx(jm))
        jlo = jm;
    else
        jhi = jm;
    end
end
return;



