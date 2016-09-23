function yp = deriv(x,y)
% function yp = deriv(x,y)
% Should be the same as idl deriv

% clearvars;
% tic;
% x = 0:0.1:10;
% y = sin(x);

n = length(x);
dx1 = x(1:end-1) - x(2:end);  % x(i) - x(i+1)
dx2 = x(1:end-2) - x(3:end);  % x(i) - x(i+2)

i = 1;
yp(i) = y(i)*(dx1(i) + dx2(i))/(dx1(i)*dx2(i)) - y(i+1)*dx2(i)/(dx1(i)*dx1(i+1)) + y(i+2)*dx1(i)/(dx2(i)*dx1(i+1));
% for i = 2:n-1
%     yp(i) = y(i-1)*dx1(i)/(dx1(i-1)*dx2(i-1)) + y(i)*(1/dx1(i) - 1/dx1(i-1)) - y(i+1)*dx1(i-1)/(dx2(i-1)*dx1(i));
% end
yp(2:n-1) = y(1:n-2).*dx1(2:n-1)./(dx1(1:n-2).*dx2(1:n-2)) + y(2:n-1)./dx1(2:n-1) - y(2:n-1)./dx1(1:n-2) - y(3:n).*dx1(1:n-2)./(dx2(1:n-2).*dx1(2:n-1));
i = n;
yp(i) = -y(i-2)*dx1(i-1)/(dx1(i-2)*dx2(i-2)) + y(i-1)*dx2(i-2)/(dx1(i-2)*dx1(i-1)) - y(i)*(dx2(i-2) + dx1(i-1))/(dx2(i-2)*dx1(i-1));
toc

% figure; hold on;
% plot(x,y,'b-.')
% plot(x,yp,'r-x')