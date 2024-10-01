function gNew = refine_gfile(g)

fac = 2;
gNew = g;

gNew.mw = (gNew.mw-1)*fac + 1;
gNew.mh = (gNew.mh-1)*fac + 1;

% Recompute r and z arrays
gNew.dR = gNew.xdim/(gNew.mw-1);
gNew.dZ = gNew.zdim/(gNew.mh-1);

gNew.r = zeros(1,gNew.mw);
gNew.z = zeros(1,gNew.mh);
for i = 0:gNew.mw-1
    gNew.r(1,i+1) = gNew.rgrid1 + gNew.dR*i;

end
for i = 0:gNew.mh-1
    gNew.z(1,i+1) = gNew.zmid - 0.5*gNew.zdim + gNew.dZ*i;
end

if isempty(gNew.cpasma)
    gNew.ip_sign = 1;
else
    gNew.ip_sign = -sign(gNew.cpasma);
end



gNew.psirz = nan(gNew.mw,gNew.mh);
for i = 1:length(gNew.z)
    gNew.psirz(:,i) = gNew.ip_sign*calc_psi(g,gNew.r,gNew.z(i)*ones(size(gNew.r)));
end

x = linspace(1,g.mw,gNew.mw);
gNew.fpol = interp1(g.fpol,x);
gNew.pres = interp1(g.pres,x);
gNew.ffprim = interp1(g.ffprim,x);
gNew.pprime = interp1(g.pprime,x);
gNew.qpsi = interp1(g.qpsi,x);

gNew = postprocess_gfile(gNew);


gNew.line1 = strrep(g.line1,strcat(' ',num2str(g.mw)),strcat(' ',num2str(gNew.mw)));

end