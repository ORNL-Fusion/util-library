function alpha_deg = calc_Bangle_g(g,R,Z)


v2 = [0,0,1];
for i = 1:length(R)
    v1 = [R(i),Z(i),0];
    vn = cross(v1,v2);
    vn = vn./norm(vn);
    b = bfield_geq_bicub(g,R(i),Z(i),0);
    b = [b.br,b.bz,b.bphi];
    b = b./norm(b);
    alpha_deg(i) = 90 - acos(dot(vn,b))*180/pi;
end




end