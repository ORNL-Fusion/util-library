function psi_bicub_coeffs = get_psi_bicub_coeffs(g)

nr = g.mw;
nz = g.mh;

psi2d = g.ip_sign*g.psirz;

dsdr = (circshift(psi2d,[-1,0])-circshift(psi2d,[1,0]))/(2*g.dR);
dsdz = (circshift(psi2d,[0,-1])-circshift(psi2d,[0,1]))/(2*g.dZ);
d2sdrdz = (circshift(psi2d,[-1,-1])-circshift(psi2d,[1,-1])-circshift(psi2d,[-1,1])+circshift(psi2d,[1,1]))/(4*g.dR*g.dZ);

bicub_mat = get_bicub_mat;
% bicub_mat_inv = get_bicub_mat_inv;
psi_bicub_coeffs = zeros(nr*nz,4,4);
for ir = 1:nr-1
    for iz = 1:nz-1
        index = iz + nz*(ir-1);
        b(:,1) = [psi2d(ir,iz),psi2d(ir+1,iz),psi2d(ir,iz+1),psi2d(ir+1,iz+1),...
             dsdr(ir,iz)*g.dR,dsdr(ir+1,iz)*g.dR,dsdr(ir,iz+1)*g.dR,dsdr(ir+1,iz+1)*g.dR,...
             dsdz(ir,iz)*g.dZ,dsdz(ir+1,iz)*g.dZ,dsdz(ir,iz+1)*g.dZ,dsdz(ir+1,iz+1)*g.dZ,...
             d2sdrdz(ir,iz)*g.dR*g.dZ,d2sdrdz(ir+1,iz)*g.dR*g.dZ,d2sdrdz(ir,iz+1)*g.dR*g.dZ,d2sdrdz(ir+1,iz+1)*g.dR*g.dZ];
         coeff = bicub_mat\b;
         psi_bicub_coeffs(index,:,1) = coeff(1:4);
         psi_bicub_coeffs(index,:,2) = coeff(5:8);
         psi_bicub_coeffs(index,:,3) = coeff(9:12);
         psi_bicub_coeffs(index,:,4) = coeff(13:16);
    end
end


% bicub_mat_inv = get_bicub_mat_inv;
% 
% psi_bicub_coeffs2 = zeros(nr*nz,4,4);
% for ir = 1:nr-1
%     for iz = 1:nz-1
%         index = iz + nz*(ir-1);
%         b(:,1) = [psi2d(ir,iz),psi2d(ir+1,iz),psi2d(ir,iz+1),psi2d(ir+1,iz+1),...
%              dsdr(ir,iz)*g.dR,dsdr(ir+1,iz)*g.dR,dsdr(ir,iz+1)*g.dR,dsdr(ir+1,iz+1)*g.dR,...
%              dsdz(ir,iz)*g.dZ,dsdz(ir+1,iz)*g.dZ,dsdz(ir,iz+1)*g.dZ,dsdz(ir+1,iz+1)*g.dZ,...
%              d2sdrdz(ir,iz)*g.dR*g.dZ,d2sdrdz(ir+1,iz)*g.dR*g.dZ,d2sdrdz(ir,iz+1)*g.dR*g.dZ,d2sdrdz(ir+1,iz+1)*g.dR*g.dZ];
%          coeff = bicub_mat_inv*b;
%          psi_bicub_coeffs2(index,:,1) = coeff(1:4);
%          psi_bicub_coeffs2(index,:,2) = coeff(5:8);
%          psi_bicub_coeffs2(index,:,3) = coeff(9:12);
%          psi_bicub_coeffs2(index,:,4) = coeff(13:16);
%     end
% end
% 
% aa=1
