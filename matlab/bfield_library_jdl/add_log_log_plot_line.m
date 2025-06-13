% clearvars;

% grab a point, then


% xx = cursor_info.Position(1);
% yy = cursor_info.Position(2);
xx = 1/200
yy = 0.2e-3
expon = 4;
x=logspace(1,3,100);
plot(x,10.^(expon*log10(x/xx)+log10(yy)),'r--','linew',1)

clear cursor_info