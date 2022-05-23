% clearvars;

% grab a point, then

expon = 6;
x=logspace(-3,2,100);
plot(x,10.^(expon*log10(x/cursor_info.Position(1))+log10(cursor_info.Position(2))),'r--','linew',1)

clear cursor_info