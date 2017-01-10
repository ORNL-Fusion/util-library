clearvars;



wout_file = 'C:\Work\Stellarator\W7X EMC3 modeling\mimic_configs\VAC\0kA\wout_w7x.nc';
wout = load_wout(wout_file);


isurf = wout.ns;


phi = 0*pi/180;

ntheta = 200;

theta = linspace(0,2*pi,ntheta);

imodes = 1:wout.mnmax;
% imodes = find(wout.xm <= 1);

for jp = 1:ntheta
    cosmn = cos(wout.xm(imodes).*theta(jp)-wout.xn(imodes)*phi);
    sinmn = sin(wout.xm(imodes).*theta(jp)-wout.xn(imodes)*phi);
    rsurf(jp) = sum(wout.rmnc(imodes,isurf).*cosmn);
    zsurf(jp) = sum(wout.zmns(imodes,isurf).*sinmn);
    if wout.asym
        rsurf(jp) = rsurf(jp) + sum(wout.rmns(imodes,isurf).*sinmn);
        zsurf(jp) = zsurf(jp) + sum(wout.zmnc(imodes,isurf).*cosmn);
    end
end


imodes = find(wout.xm <= 1 & abs(wout.xn/wout.nfp) <= 1);
for jp = 1:ntheta
    cosmn = cos(wout.xm(imodes).*theta(jp)-wout.xn(imodes)*phi);
    sinmn = sin(wout.xm(imodes).*theta(jp)-wout.xn(imodes)*phi);
    rsurf2(jp) = sum(wout.rmnc(imodes,isurf).*cosmn);
    zsurf2(jp) = sum(wout.zmns(imodes,isurf).*sinmn);
    if wout.asym
        rsurf2(jp) = rsurf2(jp) + sum(wout.rmns(imodes,isurf).*sinmn);
        zsurf2(jp) = zsurf2(jp) + sum(wout.zmnc(imodes,isurf).*cosmn);
    end
end


fprintf('Keeping %d modes\n',length(imodes))
  
for i = 1:length(imodes)         
    fprintf('RBC(%2d,%2d) = %8.4e\n',wout.xn(imodes(i))/wout.nfp,wout.xm(imodes(i)),wout.rmnc(imodes(i),isurf))
    fprintf('ZBS(%2d,%2d) = %8.4e\n',wout.xn(imodes(i))/wout.nfp,wout.xm(imodes(i)),wout.zmns(imodes(i),isurf))
end

figure; hold on; box on;
plot(rsurf.',zsurf.','r.')
plot(rsurf2,zsurf2,'b')

ves = load_W7X_vessel(0,0,phi);
plot(ves.cut.r,ves.cut.z,'k')