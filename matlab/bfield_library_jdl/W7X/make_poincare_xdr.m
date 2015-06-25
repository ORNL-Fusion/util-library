clearvars;

np = 1;
% rp = linspace(6,6.25,np);
% zp = linspace(0,0,np);
rp = 5.8;
zp = 0;
pp = 0;

ntran = 100;
distance = ntran*2*pi;

phi_period = 2*pi/5;
num_phi = distance/phi_period+1;

out_path = 'C:\Work\Stellarator\ALL_W7X_WORK\xdr_dump_read\OUTPUT\';
% fname = 'field181x181x96.w7x.1000_1000_1000_1000_+0750_+0750.vac.out';
% fname = 'fieldn_altern181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_000_s.out';
fname = 'fieldn181x181x96_w7x.0989_1010_1114_1124_+0742_-0239.09.20_022.out';
bfield_xdr = read_xdr_dump_file(out_path,fname);

% load field181x181x96.w7x.1000_1000_1000_1000_+0750_+0750.vac.bgrid.mat
% bfield_xdr = bgrid;


figure; hold on
for i = 1:np
    tic
    disp(['Running line ',num2str(i),' of ',num2str(np)])
    [r(i,:),phi(i,:),z(i,:)]=line_follow_xdr(rp(i),pp,zp(i),distance,num_phi,bfield_xdr);
%     save(['my_poincare_interp2.mat'],'r','phi','z');  
    plot(r(i,:),z(i,:),'.')
    toc
end
