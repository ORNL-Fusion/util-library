clearvars;

R1 = 1.939;
R2 = 1.872;
Z1 = -1.124;
Z2 = -1.106;


% Pappus

puffArea = 2*pi*sqrt( (R2 - R1)^2 + (Z2 - Z1)^2 )*mean([R1,R2])

