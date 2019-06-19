function graf,z,dz,pair=pair,plotm=plotm
;+
; NAME:
;      GRAF
; PURPOSE:
;      Function GRAF Returns the reduced Chi_square value
; EXPLANATION:
;        
;
; CALLING SEQUENCE:
;     GRAF, Z , DZ, PAIR=PAIR, /PLOTM
; INPUTS:
;     
;     Z - vector with the measured modulation
;     DZ - vector with the error in the modulation
;
; OPTIONAL INPUTS:
;     
;     PAIR - optional integer indentifying the source
;     PLOTM - optional integer to plot or not the modulation. If
;              PLOTM = 1 (true), the modulation and fit is ploted.
;
; OPTIONAL KEYWORD INPUTS:
;     
;
; OUTPUTS:
;     
;    The calculated Chisqrt is the output. 
;
; EXAMPLE:
;       
;  ;Chisqrt~10 (ID5617 Fild5) : badfit
;  >z =[ -0.01019  ,   0.00343  ,   0.04062  ,   0.04578   ,  0.00591  ,  -0.04084  ,  -0.05803   ,  0.03678]
;  >dz= [ 0.01161  ,   0.01134  ,   0.01143   ,  0.01132  ,  0.01147   ,  0.01168  ,   0.01133  ,   0.01207]
;  >   
;  >result = graf (z,dz) 
;  >print,result
;  > 10.262349
;
;
;  ;Chisqrt~15 (ID19972 Fild20) : goodfit
;  >z=[              -0.02779  ,  -0.03371  ,   0.02603  ,   0.02827 ,   -0.01987 ,   -0.02937 ,    0.02488  ,   0.03699]
;  >dz=[               0.00102    , 0.00100  ,   0.00100   ,  0.00117   ,  0.00100   ,  0.00101   ,  0.00100 ,    0.00100]
;  >result = graf(z,dz,plotm=1)
;  >Theoretical z:     -0.024642496    -0.032085004     0.024642499     0.032085001    -0.024642501    -0.032084999     0.024642504     0.032084997
;  >Function parameters: Q = 0.028464729, U = 0.028748205, Phi = 3.0160949
;
;  AND A WINDOW WITH THE PLOT WILL POP OUT
;
;
; PROCEDURES USED:
;       
;   CURVEFIT
;
; NOTES:
;       
; REVISON HISTORY:
;       
;-


  PAapr = dblarr(N_ELEMENTS(z)) ;Waveplate positionon
  for i=0,N_ELEMENTS(z)-1 do  PAapr[i]  = 22.5 * i 


  ;;===Fitting z===
  A = [ 0.01 , 3.14 , 0.01]     ;initial guessed parameters (Q,Phi,U)

  weigths=fltarr(N_ELEMENTS(z))
  for i=0, N_ELEMENTS(z)-1 do weigths[i]=1.
 
  zfit = CURVEFIT( PAapr, z, Weigths, A  , Sigma, FUNCTION_NAME='gfunct', CHISQ=Chisq, /DOUBLE,STATUS=status) ; [, FITA=vector] [, ITER=variable] [, ITMAX=value] [, /NODERIVATIVE] [, STATUS={0 | 1 | 2}] [, TOL=value] [, YERROR=variable] )


  ;;===uncomment to plot fit===
  if keyword_set(plotm) eq 1 then begin
     Device, RETAIN=2
     window,15,XSIZE=500,ysize=300
     plot, PAapr, z , xtitle = 'Waveplate position angle (degrees)' , ytitle = 'Amplitud of modulation (z)', PSYM = 2 , xrange=[-10,200],background=0,COLOR=16777215
     oploterror , PAapr , z , dz
     X=FLOAT(INDGEN(180))*!pi/180
     oplot, X*180/!pi , A[0] * cos (4*X  + A[1])  + A[2] * sin ( 4*X  + A[1]) , LINESTYLE=2
     ztheo =  A[0] * cos (4*PAapr*!pi/180  + A[1])  + A[2] * sin ( 4*PAapr*!pi/180  + A[1])
     print,'Theoretical z: ', ztheo
     print,'Function parameters: Q = ',STRTRIM(A[0],2),', U = ', STRTRIM(A[2],2),', Phi = ',STRTRIM(A[1],2)
  endif

  ;;===CALUCLATING reduced CHISQRT===
  Df =  N_ELEMENTS(z) - N_ELEMENTS(A)    ;degrees of freedom= Ndata - Nfitting_parameters
  Chisqrt = 1./Df*total((z-zfit)^2/dz^2) ;reduced Chisqrt
  ;Chisqrt_crit = 1./Df*(1.^2) * N_ELEMENTS(z) ;Critical reduced chisqrt value when: (z-zfit)=1.dz 
 

  ;Chisq =  1./Df*total((z-zfit)^2) ;chisq from Curvefit.

  ;;===printing (or not) the pair number when fail to converge===
  if keyword_set(pair) eq 1 then begin
     if status ne 0 then print,'Pair whose modulation failed to converge: ',pair
  endif else begin
     if status ne 0 then print,'Modulation failed to converge.'
  endelse

     
  return,Chisqrt
  
end
