function ves_cut = cut_W7X_vessel(ves,thisphi_rad)
% if thisphi_rad [rad] is provided then a single cut is returned interpolated
% at this angle

if thisphi_rad == 2*pi/ves.nsym    
    thisphi_fp1 = thisphi_rad;
else
    thisphi_fp1 = mod(thisphi_rad,2*pi/ves.nsym);
end


ves_cut.x = zeros(1,ves.npol);
ves_cut.y = zeros(1,ves.npol);
ves_cut.z = zeros(1,ves.npol);
ves_cut.r = zeros(1,ves.npol);
ves_cut.phi = zeros(1,ves.npol);
for j=1:ves.npol
    ves_cut.z(j) = interp1(ves.phi(:,j),ves.z(:,j),thisphi_fp1);
    ves_cut.r(j) = interp1(ves.phi(:,j),ves.r(:,j),thisphi_fp1);
end
ves_cut.phi(:) = thisphi_rad;
ves_cut.x = ves_cut.r.*cos(thisphi_rad);
ves_cut.y = ves_cut.r.*sin(thisphi_rad);


