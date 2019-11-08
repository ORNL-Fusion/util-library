function y = evaluate_tanh_fit(c,x,param)
%
% ; tanh function with cubic or quartic inner and linear
% ;                                to quadratic outer extensions and derivative=0 at param
% ;       0.5*(c[2]-c[3])* ( pz1*exp(z) - pz2*exp(-z) )/
% ;                                     ( exp(z) + exp(-z) ) + 0.5*(c[2]+c[3])
% ;          where z = 2.*(c[0]-x)/c[1]
% ;      if param = None:
% ;          pz1 = 1.+ c[4]*z + c[5]*z*z + c[6]*z*z*z
% ;       else:
% ;           pz1 = 1 + cder*z + c[4]*z*z + c[5]*z*z*z + c[6]*z*z*z*z
% ;           where cder = -( 2.0*c[4]*z0 + + 3.0*c[5]*z0*z0 + 4.0*c[6]*z0*z0*z0 )
% ;           and z0 = 2.*( c[0] - param )/c[1]
% ;       pz2 = 1 + ( c[7]*z + (c[8]*z*z) ) depending on whether there are 7,8,or 9
% ;             coefficients specified
% ;       c0 = SYMMETRY POINT
% ;       c1 = FULL WIDTH
% ;       c2 = HEIGHT
% ;       c3 = OFFSET
% ;       c4 = SLOPE OR QUADRATIC (IF ZERO DER) INNER
% ;       c5 = QUADRADIC OR CUBIC (IF ZERO DER) INNER
% ;       c6 = CUBIC OR QUARTIC (IF ZERO DER) INNER
% ;       c7 = SLOPE OUTER
% ;       c8 = QUADRATIC OUTER

if nargin < 3
    param = [];
else
    if param ~= 0 && ~isempty(param)
        param = 1;
    end
end

c0 = c(1); c1 = c(2); c2 = c(3); c3 = c(4); c4 = c(5); 
if length(c) > 5
    c5 = c(6); 
end
if length(c) > 6 
    c6 = c(7); 
end
if length(c) > 7
    c7 = c(8); 
else
    c7 = 0;
end
if length(c) > 8
    c8 = c(9);
else 
    c8 = 0;
end

z = 2*(c0-x)/c1;
ez = exp(z);
emz = exp(-z);
if length(c) == 5
    y = 0.5*(c2-c3)*((1+c4*z).*ez - emz)./(ez + emz) + 0.5*(c2+c3);
    return;
end

if length(c) == 6
    if ~isempty(param)
        z0 = 2*( c0 - param )/c1;
        cder = -(2*c3*z0 + 3*c4*z0^2 + 4*c5*z0^2);
        pz1 = 1 + cder*z + c3*z.^2 + c4*z.^3 + c5*z.^4;
    else
        pz1 = 1 + c3*z + c4*z.^2 + c5*z.^3;        
    end
    y = 0.5*c2*(pz1.*ez - emz)./(ez + emz) + 0.5*c2;
    return;
end

if ~isempty(param)
    z0 = 2*( c0 - param )/c1;
    cder = -(2*c4*z0 + 3*c5*z0^2 + 4*c6*z0^2);
    pz1 = 1 + cder*z + c4*z.^2 + c5*z.^3 + c5*z.^4;
else
    pz1 = 1 + c4*z + c5*z.^2 + c6*z.^3;
    pz2 = 1 + c7*z + c8*z.^2;
end
y = 0.5*(c2-c3)*(pz1.*ez - pz2.*emz)./(ez + emz) + 0.5*(c2+c3);

