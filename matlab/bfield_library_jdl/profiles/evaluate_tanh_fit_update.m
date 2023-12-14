function y = evaluate_tanh_fit_update(coeffs,x,type,param)
% 
% type can be 
%  'hy1a'
%  'hy1b'
%  'tanh'
%  'tanh0'
%
% param only used for tanh type

if nargin < 4
    param = [];
end


switch lower(type)
    case 'hy1a'
        % HY1A
        %     ---------------------------------------------------------------------
        %     Piecewise linear function forming a single hyperbola
        %       Region 1, x > ( c[3] - c[1] ) / c[2], y = c[3]
        %       Region 2, x < ( c[3] - c[1] ) / c[2], y = c[1] + c[2] * x = v
        %
        %       v = c[1] + c[2] * x
        %
        %       y = 0.5 * ( v  + numpy.sqrt( v**2 + 4.0 * c[0]**2 ) ) + c[3]
        %
        %       c[0] = curvature parameter - corner of hyperbola joining two linear segments
        %       c[1] = y-intercept of line in region 2
        %       c[2] = slope of line in region 2
        %       c[3] = offset in region 1
        %     ---------------------------------------------------------------------

        c = coeffs;
        v = c(2) + c(3) * x;
        y = 0.5 * ( v + sqrt( v.^2 + 4.0 * c(1)^2 ) ) + c(4);

    case 'hy1b'
        % HY1B
        %     ---------------------------------------------------------------------
        %     Piecewise linear function forming a single hyperbola
        %       Region 1, - ( c[1] - c[3] ) / ( c[2] - c[4] ) < x, y = c[1] + c[2] * x = u
        %       Region 2, x < - ( c[1] - c[3] ) / ( c[2] - c[4] ), y = c[3] + c[4] * x = v
        %
        %       u = c[1] + c[2] * x
        %       v = c[3] + c[4] * x
        %
        %       y = 0.5 * ( ( u + v ) + numpy.sqrt( ( u - v )**2 + 4.0 * c[0]**2 ) )
        %
        %       c[0] = curvature parameter - corner of hyperbola joining two linear segments
        %       c[1] = y-intercept of line in region 1
        %       c[2] = slope of line in region 1
        %       c[3] = y-intercept of line in region 2
        %       c[4] = slope of line in region 2
        %     ---------------------------------------------------------------------

        c = coeffs;
        u = c(2) + c(3) * x;
        v = c(4) + c(5) * x;

        y = 0.5 * ( (u + v) + sqrt( (u - v).^2 + 4.0 * c(1)^2 ) );

    case 'tanh'
        % TANH
        %
        % ---------------------------------------------------------------------
        % tanh_multi(c,x,param=None): tanh function with cubic or quartic inner and linear
        %                             to quadratic outer extensions and derivative=0 at param
        %    0.5*(c[2]-c[3])* ( pz1*exp(z) - pz2*exp(-z) )/
        %                                  ( exp(z) + exp(-z) ) + 0.5*(c[2]+c[3])
        %       where z = 2.*(c[0]-x)/c[1]
        %   if param = None:
        %       pz1 = 1.+ c[4]*z + c[5]*z*z + c[6]*z*z*z
        %   else:
        %       pz1 = 1 + cder*z + c[4]*z*z + c[5]*z*z*z + c[6]*z*z*z*z
        %       where cder = -( 2.0*c[4]*z0 + + 3.0*c[5]*z0*z0 + 4.0*c[6]*z0*z0*z0 )
        %       and z0 = 2.*( c[0] - param )/c[1]
        %   pz2 = 1 + ( c[7]*z + (c[8]*z*z) ) depending on whether there are 7,8,or 9
        %         coefficients specified
        %   c0 = SYMMETRY POINT
        %   c1 = FULL WIDTH
        %   c2 = HEIGHT
        %   c3 = OFFSET
        %   c4 = SLOPE OR QUADRATIC (IF ZERO DER) INNER
        %   c5 = QUADRADIC OR CUBIC (IF ZERO DER) INNER
        %   c6 = CUBIC OR QUARTIC (IF ZERO DER) INNER
        %   c7 = SLOPE OUTER
        %   c8 = QUADRATIC OUTER
        %
        % ---------------------------------------------------------------------
        
        c = zeros(1,9);
        c(1:length(coeffs)) = coeffs;

        if isempty(param)
            z = 2 * (c(1) - x) / c(2);
            pz1 = 1 + c(5)*z + c(6)*z.^2 + c(7)*z.^3;
            pz2 = 1 + (c(8)*z + c(9)*z.^2);
        else
            z0 = 2 * (c(1) - param) / c(2);
            cder = -(2.0 * c(5) * z0 + 3.0 * c(6) * z0.^2 + 4.0 * c(7) * z0.^3);
            z = 2 * (c(1) - x) / c(2);
            pz1 = 1 + cder*z + c(5)*z.^2 + c(6)*z.^3 + c(7)*z.^4;
            pz2 = 1 + (c(8)*z + (c(9)*z.^2));
        end
        y = 0.5 * (c(3) - c(4)) * (pz1 .* exp(z) - pz2 .* exp(-z)) ./ (exp(z) + exp(-z)) + 0.5 * (c(3) + c(4));     

    case 'tanh0'
        % TANH0
        %
        %
        % -----------------------------------------------------------------------
        % tnh0(c,x,param=None): new tanh function
        %    0.5*(c[2]-c[3])* ( (1+c[4]*z)*exp(z) - exp(-z) )/( exp(z) + exp(-z) ) + 0.5*(c[2]+c[3])
        %       where z = 2.*(c[0]-x)/c[1]
        %    c[0] = symmetry point
        %    c[1] = full width
        %    c[2] = height
        %    c[3] = offset
        %    c[4] = alpha
        % ----------------------------------------------------------------------

        c = coeffs;
        z = 2 * (c(1) - x) / c(2);
        y = 0.5 * (c(3) - c(4)) * ( (1 + c(5) * z) .* exp(z) - exp(-z) ) ./ ( exp(z) + exp(-z) ) + 0.5 * (c(3) + c(4));
        

    otherwise
        error('Bad type')
end

