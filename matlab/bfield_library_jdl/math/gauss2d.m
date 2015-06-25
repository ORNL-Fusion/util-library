function f=gauss2d(x,y,A,x0,y0,sigx,sigy)
% function f=gauss2d(x,y,A,x0,y0,sigx,sigy)
% Returns evaluation of 2D gaussian
% f(x,y) = A*exp( -( (x-x0)^2/(2*sigx^2) +(y-y0)^2/(2*sigy^2)  ))
% V = 2*pi*sigx*sigy  (total volume under curve)

% for i = 1:length(x)
%     for j = 1:length(y)
        f = A*exp( -( (x-x0).^2./(2*sigx^2) +(y-y0).^2./(2*sigy^2)  ));
%     end
% end