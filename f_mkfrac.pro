; f_mkfrac(x): 
; 
; A function for converting an integral or half-integral number to 
; a string with the corresponding fractional number.
;   
;   inp: any integral or half-integral number x with abs(x) < 50
;   out: 5 char right adjusted string containing a fractional representation of x, 
;           e.g. " -3/2", "-11/2" or "    2"
;   
;   Jon Grumer, Uppsala, 2018

function f_mkfrac, x

  ; Error handling: check if x is within range for printing in a 5 char string
  if abs(x) ge 50 then message, 'Input value too large, abs(x) has to be <= 50'

  ; Check if x = 0
  if x eq 0. then begin
    s = '    0'
  ; Check if x is integral  
  endif else if abs(x*2) mod 2 eq 0 then begin
    s = string(x, '(i5)')
  ; Check if x is half-integral
  endif else begin
    ; Check if 2x is larger than 10
    if abs(x*2) ge 10 then begin
      s1 = string(x*2, '(i3)')
    endif else begin
      s1 = " " + string(x*2, '(i2)')
    endelse
     s = s1 + '/2'
  endelse

  return, s

end