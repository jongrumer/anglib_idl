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
;   Dependencies: f_wig3j.pro
;
;   Jon Grumer, Uppsala University, 2018

function f_cg, j1, j2, j3, m1, m2, m3
  compile_opt idl2

  ; Determine the Clebsch-Gordan coeff via Wigner 3j symbol.

  cg = (-1.0d)^(j1-j2+m3) * sqrt(2.0d*double(j3) + 1.0d) * f_wig3j(j1,j2,j3,m1,m2,-m3)

  return, cg

end
