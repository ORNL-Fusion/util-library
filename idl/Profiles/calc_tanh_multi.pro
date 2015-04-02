function calc_tanh_multi,c,x,param=param

; tanh function with cubic or quartic inner and linear
;                                to quadratic outer extensions and derivative=0 at param
;       0.5*(c[2]-c[3])* ( pz1*exp(z) - pz2*exp(-z) )/
;                                     ( exp(z) + exp(-z) ) + 0.5*(c[2]+c[3])
;          where z = 2.*(c[0]-x)/c[1]
;      if param = None:
;          pz1 = 1.+ c[4]*z + c[5]*z*z + c[6]*z*z*z
;       else:
;           pz1 = 1 + cder*z + c[4]*z*z + c[5]*z*z*z + c[6]*z*z*z*z
;           where cder = -( 2.0*c[4]*z0 + + 3.0*c[5]*z0*z0 + 4.0*c[6]*z0*z0*z0 )
;           and z0 = 2.*( c[0] - param )/c[1]
;       pz2 = 1 + ( c[7]*z + (c[8]*z*z) ) depending on whether there are 7,8,or 9
;             coefficients specified
;       c0 = SYMMETRY POINT
;       c1 = FULL WIDTH
;       c2 = HEIGHT
;       c3 = OFFSET
;       c4 = SLOPE OR QUADRATIC (IF ZERO DER) INNER
;       c5 = QUADRADIC OR CUBIC (IF ZERO DER) INNER
;       c6 = CUBIC OR QUARTIC (IF ZERO DER) INNER
;       c7 = SLOPE OUTER
;       c8 = QUADRATIC OUTER

z = 2.*(c[0]-x)/c[1]
if n_elements(c) eq 5 then begin
    out= 0.5*(c[2]-c[3])* ( (1+c[4]*z)*exp(z) - exp(-z) )/( exp(z) + exp(-z) ) + 0.5*(c[2]+c[3])
endif else begin

    if n_elements(c) eq 6 then begin
        if keyword_set(param) then begin
            z0 = 2.*( c[0] - param )/c[1]
            cder = -( 2.0*c[3]*z0 + 3.0*c[4]*z0*z0 + 4.0*c[5]*z0*z0*z0 )
            pz1 = 1 + cder*z + c[3]*z*z + c[4]*z*z*z + c[5]*z*z*z*z
        endif else pz1 = 1.+ c[3]*z + c[4]*z*z + c[5]*z*z*z
        out =  0.5*c[2]* ( pz1*exp(z) - exp(-z) )/( exp(z) + exp(-z) ) + 0.5*c[2]
    endif else begin

        if keyword_set(param) then begin
            z0 = 2.*( c[0] - param )/c[1]
            cder = -( 2.0*c[4]*z0 + 3.0*c[5]*z0*z0 + 4.0*c[6]*z0*z0*z0 )
            pz1 = 1 + cder*z + c[4]*z*z + c[5]*z*z*z + c[6]*z*z*z*z
        endif else pz1 = 1.+ c[4]*z + c[5]*z*z + c[6]*z*z*z
        pz2 = 1
        if n_elements(c) gt 7 then pz2 += c(7)*z
        if n_elements(c) gt 8 then pz2 += c(8)*z*z

        out =  0.5*(c[2]-c[3])* ( pz1*exp(z) - pz2*exp(-z) )/ $
          ( exp(z) + exp(-z) ) + 0.5*(c[2]+c[3])
    endelse
endelse

return,out
end
