function styles = styflipper(n)

s = {'-','--',':','-.'};
ns = length(s);

for i = 1:n
    styles{i} = s(mod(i-1,ns)+1);
end