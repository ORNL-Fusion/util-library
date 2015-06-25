function [ current ] = set_w7x_current( numelc,taper )

% taper = [I1,I2,I3,I4,I5,IA,IB];
current = zeros(sum(numelc),1);

t2 = [taper(1:5),taper(5:-1:1),taper([6,7,7,6])];
t2 = [t2,t2,t2,t2,t2];

icoil = 1;
icount = 1;
for i = 1:sum(numelc)
%     if i == 971
%         a=1;
%     end
    if icount < numelc(icoil)
        current(i) = t2(icoil);
        icount = icount + 1;
    elseif icount == numelc(icoil)
        current(i) = 0.d0;
        icount = 1;
        icoil = icoil+1;
    end
    
end

% 
% for i=1:length(taper)
%     n = numelc(i);
%     current(1:n) = taper(i);
% end

end
