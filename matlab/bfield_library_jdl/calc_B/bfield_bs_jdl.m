function [Bx,By,Bz,Btot]=bfield_bs_jdl(P_x,P_y,P_z,coil,current)
% [Bx,By,Bz,Btot]=bfield_bs_jdl(P_x,P_y,P_z,coil,current)

I=current.*1e-7; %this is mu0*I/4pi
npts = length(P_x);
Bx=zeros(npts,1);
By=zeros(npts,1);
Bz=zeros(npts,1);

coils_x = coil(:,1);
coils_y = coil(:,2);
coils_z = coil(:,3);
nfil=length(coils_x)-1;

for i = 1:npts

    Ri_x = P_x(i) - coils_x(1:nfil);     %R_i is from beg of filament to point x
    Ri_y = P_y(i) - coils_y(1:nfil);
    Ri_z = P_z(i) - coils_z(1:nfil);
    Ri_mag = sqrt(Ri_x.*Ri_x + Ri_y.*Ri_y + Ri_z.*Ri_z);
    
    Rf_x=[Ri_x(2:nfil); P_x(i) - coils_x(nfil+1)];                    %R_f is from end of filament to x
    Rf_y=[Ri_y(2:nfil); P_y(i) - coils_y(nfil+1)];                    %it is thus the same as Ri except one point before
    Rf_z=[Ri_z(2:nfil); P_z(i) - coils_z(nfil+1)];
    Rf_mag=[Ri_mag(2:nfil);sqrt(Rf_x(nfil)*Rf_x(nfil) + Rf_y(nfil)*Rf_y(nfil) + Rf_z(nfil)*Rf_z(nfil))];
    
    RicrRf_x=Ri_y.*Rf_z-Ri_z.*Rf_y;      %cross product: RixRf
    RicrRf_y=Ri_z.*Rf_x-Ri_x.*Rf_z;
    RicrRf_z=Ri_x.*Rf_y-Ri_y.*Rf_x;
    
    RidotRf = Ri_x.*Rf_x + Ri_y.*Rf_y + Ri_z.*Rf_z;
    
    back = I(1:nfil).*(Ri_mag+Rf_mag)./(Ri_mag.*Rf_mag.*(Ri_mag.*Rf_mag+RidotRf));
    
    Bx(i)=sum(RicrRf_x.*back);
    By(i)=sum(RicrRf_y.*back);
    Bz(i)=sum(RicrRf_z.*back);       
    
end

if nargin > 3
    Btot = sqrt(Bx.^2+By.^2+Bz.^2);
end

end
