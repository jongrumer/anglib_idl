;
;   IDL Function for computing Wigner 6j-symbols: {j1 j2 j3}
;   ============================================  {l1 l2 l3}
;
;   Based on and follows the notation of expression (5.23) in
;   R.D. Cowan, "The Theory of Atomic Structure and Spectra"
;
;   Jon Grumer, Uppsala University, 2018

function f_wig6j, j1, j2, j3, l1, l2, l3
  compile_opt idl2

  if ( j1 ge 0 and j2 ge 0 and j3 ge 0 and l1 ge 0 and l2 ge 0 and l3 ge 0 and $ ; all ji and li >= 0
    4*j1 mod 2 eq 0 and 4*j2 mod 2 eq 0 and 4*j3 mod 2 eq 0 and $                ; all ji and li must be integral...
    4*l1 mod 2 eq 0 and 4*l2 mod 2 eq 0 and 4*l3 mod 2 eq 0 and $                ;               or half integral
    2*(j1+j2+j3) mod 2 eq 0 and $                                                ; j1+j2+j3 integral    
    2*(j1+l2+l3) mod 2 eq 0 and $                                                ; j1+l2+l3 integral    
    2*(l1+j2+l3) mod 2 eq 0 and $                                                ; l1+j2+l3 integral
    2*(l1+l2+j3) mod 2 eq 0 and $                                                ; l1+l2+j3 integral
    j1 + j2 ge j3 and j2 + j3 ge j1 and j3 + j1 ge j2 and $                      ; (j1,j2,j3) triangle relations (E.g. Cowan (5.26))
    j1 + l2 ge l3 and l2 + l3 ge j1 and l3 + j1 ge l2 and $                      ; (j1,l2,l3) triangle relations  
    l1 + j2 ge l3 and j2 + l3 ge l1 and l3 + l1 ge j2 and $                      ; (j1,j2,l3) triangle relations
    l1 + l2 ge j3 and l2 + j3 ge l1 and j3 + l1 ge l2 $                          ; (l1,l2,j3) triangle relations
    ) then begin

    kmin = max([j1+j2+j3, j1+l2+l3, l1+j2+l3, l1+l2+j3])
    kmax = min([j1+j2+l1+l2, j2+j3+l2+l3, j3+j1+l3+l1])

    w6j = 0.0d
    
    for k = kmin, kmax do begin
      numerator   = factorial(k+1)
      denominator = factorial(k-j1-j2-j3) * factorial(k-j1-l2-l3) * $
                    factorial(k-l1-j2-l3) * factorial(k-l1-l2-j3) * $
                    factorial(j1+j2+l1+l2-k) * factorial(j2+j3+l2+l3-k) * factorial(j3+j1+l3+l1-k)
      w6j = w6j + (-1.0d)^k * (double(numerator) / double(denominator))
    end
 
    w6j = w6j * f_fac(j1,j2,j3) * f_fac(j1,l2,l3) * f_fac(l1,j2,l3) * f_fac(l1,l2,j3)
    
  endif else begin
    w6j = 0.0d
  end
  
  return, w6j

end

function f_fac, a, b, c
  compile_opt idl2, hidden
  return, sqrt( factorial(a+b-c) * factorial(a-b+c) *  $ 
                factorial(-a+b+c) / factorial(a+b+c+1) )
end