function syms = symflipper(n)

s = {'o','*','s','<','d','p','h','v'};
ns = length(s);

for i = 1:n
    syms{i} = char(s(mod(i-1,ns)+1));
end