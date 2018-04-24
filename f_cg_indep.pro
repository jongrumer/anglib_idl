;
;   IDL Function for computing Clebsch-Gordan Coefficients: [j1 j2 j3]
;   ======================================================  [m1 m2 m3]
;
;   Based on and follows the notation of expression (5.1) in 
;   R.D. Cowan, "The Theory of Atomic Structure and Spectra"
;   which gives the Wigner 3j-symbol (W3J). The CG coefficient
;   is then determined from

;   CG(j1,j2,j3,m1,m2,m3) = (-1)^(j1-j2+m3) * sqrt(2*j3 + 1) ...
;                           * W3J(j1,j2,j3,m1,m2,-m3)
;   
;   Where the sign change in m3 should be noted.
;   
;   See e.g. section 8.1.2 of "Quantum Theory of Angular Momentum" 
;   by Varshalovich, Moskalev and Khersonskii.
;
;   This version is fully independent of any external routines.
;
;   Jon Grumer, Uppsala University, 2018

function f_cg_indep, j1, j2, j3, m1, m2, m3_tmp
  compile_opt idl2
  
  ; Flip the sign of m3, to determine the correct 3j-symbol.
  m3 = - m3_tmp

  if ( j1 ge abs(m1)        and j2 ge abs(m2)           and j3 ge abs(m3)           and $ ; ji >= abs(mi) >= 0, for all i's
    4*j1 mod 2 eq 0         and 4*j1 mod 2 eq 0         and 4*j1 mod 2 eq 0         and $ ; ji integral or half-integral, for all i's
    2*(j1-m1) mod 2 eq 0    and 2*(j2-m2) mod 2 eq 0    and 2*(j3-m3) mod 2 eq 0    and $ ; both ji and mi must ebe either integral, for all i's
    2*(j1+m1) mod 2 eq 0    and 2*(j2+m2) mod 2 eq 0    and 2*(j3+m3) mod 2 eq 0    and $ ; or half-integral
    2*(m1+m2+m3) mod 2 eq 0 and 2*(j1-j2-m3) mod 2 eq 0 and 2*(j1+j2+j3) mod 2 eq 0 and $ ; m1+m2+m3, j1-j2-m3 and j1+j2+j3 == integral
    2*(j1+j2) ge 2*j3       and 2*(j2+j3) ge 2*j1       and 2*(j3+j1) ge 2*j2       and $ ; (j1,j2,j3) triangle relations
    abs(2*(m1+m2+m3)) eq 0 $
    ) then begin

    kmin = max( [0, j2-j3-m1, j1-j3+m2])
    kmax = min( [j1+j2-j3, j1-m1, j2+m2])
    
    w3j = 0.0d
    
    for k = kmin, kmax do begin
      numerator = 1.0d
      denominator = factorial(k) * factorial(j1+j2-j3-k) * factorial(j1-m1-k) *        $
                    factorial(j2+m2-k) * factorial(j3-j2+m1+k) * factorial(j3-j1-m2+k)
      w3j = w3j + (-1.0d)^k * double(numerator) / double(denominator)
    end
    numerator = factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) *     $
                factorial(j1-m1) * factorial(j1+m1) * factorial(j2-m2) *               $
                factorial(j2+m2) * factorial(j3-m3) * factorial(j3+m3)
    denominator = factorial(j1+j2+j3+1)
    w3j = w3j * (-1.0d)^(j1-j2-m3) * sqrt( double(numerator) / double(denominator) )

    ; Determine the Clebsch-Gordan coeff from the calculated Wigner 3j symbol with flipped sign of m3.
    cg  = w3j * (-1.0d)^(j1-j2+m3_tmp) * sqrt(2.0d*j3 + 1.0d)
  
  endif else begin
    cg = 0.0d
  end

  return, cg

end
