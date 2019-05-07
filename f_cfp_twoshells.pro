;
;   F_CFP_TWOSHELLS
;
;   Jon Grumer, Uppsala University, 2019
;
;   Function for determination of mixed, two-shell CFP's used to build
;   antisymmetric atomic states from states with two involved subshells
;   with various cores and parent terms:
;
;   Psi( [l_p^w_p L_p S_p, l_a^w_a L_a S_a] LS)
;
;     = CFP x psi( [ (l_p^w_p L_p S_p,  l_a^w_a L_ac S_ac) L_c S_c ] l_as     LS )
;                    passive            active             core      active   total
;                    shell              shell              term      electron term
;
;   Parameters (in input order)
;   ===========================
;   cfp_s         single-shell cfp of the active shell (l_a^w_a-1 |} l_a^w_a)
;   L_p,  S_p       passive shell term
;   L_ac, S_ac      active shell core (parent) term
;   L_a,  S_a       active shell term
;   L_c,  S_c       core term
;   l_as            active single electron l (and s = 0.5)
;   L,    S         total term
;
;   Dependencies:
;   =============
;   - f_cg_indep.pro
;
;   Test cases (copy-paste the function calls):
;   ===========================================
;   Case 1: [ (3d2 3F, 4s 2S) 4F ] 4s 3F --> 0.8164965 or sqrt(4/6)
;                      cfp_s  L_p S_p   L_ac S_ac   L_a S_a   L_c S_c   l_as  L  S
;   >> f_cfp_twoshells(1,       3,  1,    0,   0.5,   0,  0,    3,  1.5,  0,    3, 1)

;   Case 2: [ (3d2 3F, 4s 2S) 2F ] 4s 3F --> -0.577350 or -sqrt(2/6)
;                      cfp_s  L_p S_p   L_ac S_ac   L_a S_a   L_c S_c   l_as  L  S
;   >> f_cfp_twoshells(1,       3,  1,    0,   0.5,   0,  0,    3,  0.5,  0,    3, 1)
;

function f_cfp_twoshells, cfp_s, L_p, S_p,  L_ac, S_ac,  L_a, S_a,  L_c, S_c,  l_as,  L, S
   compile_opt idl2

   ; normalization factor (2L+1)(2S+1)
   norm_fac = ((2d0*double(L)+1d0)*(2d0*double(S)+1d0))

   ; single-electron spin
   s_as = 0.5d0

   ; for a given active shell core term (L_ac, S_ac) the two shell-cfp is:

   ;print
   ;print, '   Clebsh-Gordan coefficients'
   ;print, '   L1         S1         L_a         S_a         L_c         S_c         L2         S2         SUM'

   sum = 0d0
   for ML = -L, L do begin
      for MS = -S, S do begin

         for ML_p = -L_p, L_p do begin
            for MS_p = -S_p, S_p do begin
               for ML_a = -L_a, L_a do begin
                  for MS_a = -S_a, S_a do begin

                     ; decoupling active and passive subshell
                     CG_L1 = f_cg_indep(L_a,L_p,L, ML_a,ML_p,ML)
                     CG_S1 = f_cg_indep(S_a,S_p,S, MS_a,MS_p,MS)

                     one   = CG_L1*CG_S1

                     for ML_ac = -L_ac, L_ac do begin
                        for MS_ac = -S_ac, S_ac do begin
                           for ml_as = -l_as, l_as do begin
                              for ms_as = -s_as, s_as do begin

                                 ; decoupling active electron
                                 CG_L_a = f_cg_indep(L_ac,l_as,L_a, ML_ac,ml_as,ML_a)
                                 CG_S_a = f_cg_indep(S_ac,s_as,S_a, MS_ac,ms_as,MS_a)

                                 two = CG_L_a*CG_S_a

                                 for ML_c = -L_c, L_c do begin
                                    for MS_c = -S_c, S_c do begin

                                       ;  re-couple passive shell to remaining electrons in
                                       ;  active shell to a mixed shell core
                                       CG_L_c = f_cg_indep(L_p,L_ac,L_c, ML_p,ML_ac,ML_c)
                                       CG_S_c = f_cg_indep(S_p,S_ac,S_c, MS_p,MS_ac,MS_c)

                                       three = CG_L_c*CG_S_c

                                       ; couple new core to decoupled active electron
                                       CG_L2 = f_cg_indep(L_c,l_as,L, ML_c,ml_as,ML)
                                       CG_S2 = f_cg_indep(S_c,s_as,S, MS_c,ms_as,MS)

                                       four = CG_L2*CG_S2

                                       ; if all three contributions are non-zero, then add to cfp sum
                                       if abs(one*two*three*four) gt 0d0 then begin
                                          sum = sum + one*two*three*four

                                          ;print, CG_L1, CG_S1, CG_L_a, CG_S_a, CG_L_c, $
                                          ;  CG_S_c,CG_L2, CG_S2, cfp, $
                                          ;  format='(9f11.6)'

                                       endif

                                    endfor
                                 endfor
                              endfor
                           endfor
                        endfor
                     endfor
                  endfor
               endfor
            endfor
         endfor
      endfor
   endfor

   ; multiply with single-shell cfp and normalize
   cfp = sum * cfp_s / norm_fac

   ; print two-shell cfp
   ;
   ;terms=['S','P','D','F','G','H','I']
   ;print
   ;print, '   CFP [(',(2*S_p+1),terms[fix(L_p)],',',(2*S_ac+1),terms[fix(L_ac)],')', $
   ;   (2*S_c+1),terms[fix(L_c)],']',(2*S+1),terms[fix(L)],' = ', $
   ;   cfp, $
   ;   format='(a9,i1,2a1,i1,2a1,i1,2a1,i1,a1,a3,f9.6)'
   ;print

   return, cfp

end
