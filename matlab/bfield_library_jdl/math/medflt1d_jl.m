function xout = medflt1d_jl(x,m,pad,omitnan)
% function xout = medflt1d_jl(x,m,pad,omitnan)
% x is 1d array, m is number of points over which to apply median
% end points are padded (or make pad = 'zero')

if nargin < 3
    pad = 'end';
end
if nargin < 4
    omitnan = 0;
end

% % % m = 10;
% % % pad = 'zero';
% % % fs = 100;
% % % t = 0:1/fs:1;
% % % x = sin(2*pi*t*3)+0.25*sin(2*pi*t*40);
% % % 

if omitnan == 1
    miss = 'omitnan';
elseif omitnan == 0
    miss = 'includenan';
end

iflip = 0;
if iscolumn(x) 
    x = x.';
    iflip = 1;
end


if ~mod(m,2);
    dimin = m/2;
    dimax = m/2+1;    
else
    dimin = (m-1)/2;
    dimax = dimin;
end

x2 = zeros(1,length(x)+dimin+dimax);
x2(dimin+1:end-dimax) = x;
if strcmp(pad,'end')        
    x2(1:dimin) = x(1);
    x2(end-dimax+1:end) = x(end);
elseif strcmp(pad,'zero')
    x2(1:dimin) = 0;
    x2(end-dimax+1:end) = 0;
end    
xout = x2;

for i = dimin+1:length(x2)-dimax
    xout(i) = median(x2(i-dimin:i+dimax),2,miss);
end

xout= xout(dimin+1:end-dimax);

% % figure; hold on; box on;
% % plot(t,x,'b')
% % plot(t,xout,'r')

if iflip
    xout = xout.';
end