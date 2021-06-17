function alpha_deg = calc_Bangle_g(g,R1,Z1,R2,Z2)
%  alpha_deg = calc_Bangle_g(g,R1,Z1,R2,Z2)

v2 = [0,0,1];
for i = 1:length(R1)
    v1 = [R2(i)-R1(i),Z2(i)-Z1(i),0];
    vn = cross(v1,v2);
    vn = vn./norm(vn);
    b = bfield_geq_bicub(g,0.5*(R2(i)+R1(i)),0.5*(Z2(i)+Z1(i)),0);
    b = [b.br,b.bz,b.bphi];
    b = b./norm(b);
    alpha_deg(i) = 90 - acos(dot(vn,b))*180/pi;
end




end