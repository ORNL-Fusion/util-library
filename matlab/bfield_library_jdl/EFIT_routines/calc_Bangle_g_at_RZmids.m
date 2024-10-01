function Bangle = calc_Bangle_g_at_RZmids(g,R,Z)
% calc_Bangle_g_at_RZmids(g,R,Z)
% Just a wrapper to pass R,Z arrays to calc_Bangle_g
% output will be one smaller because Bangle is at midpoints
% JDL 2024

for i = 1:length(R) - 1
    RZ1 = [R(i),Z(i)];
    RZ2 = [R(i+1),Z(i+1)];
    BangleTmp = calc_Bangle_g(g,RZ1,RZ2,1);
    Bangle.alpha_deg(i) = BangleTmp.alpha_deg;
    Bangle.beta_deg(i) = BangleTmp.beta_deg;
end
end
