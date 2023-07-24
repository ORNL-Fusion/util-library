function area = pappus(R1,Z1,R2,Z2)

% R1 = 2
% R2 = 1.9;
% Z1 = -1;
% Z2 = -1.;


% Pappus

area = 2*pi*sqrt( (R2 - R1)^2 + (Z2 - Z1)^2 )*mean([R1,R2])

