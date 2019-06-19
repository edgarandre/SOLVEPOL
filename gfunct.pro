;;function to fit
PRO gfunct, PAapr, A, F , pder

  X = PAapr *!pi/180.

  F = A[0] * cos ( 4*X + A[1]) + A[2] * sin ( 4*X + A[1]) 

  dA0 = cos (4*X + A[1])
  dA1 = - A[0] * sin (4*X + A[1]) + A[2] * cos ( 4*X + A[1]) 
  dA2 = sin (4*X + A[1]) 


  pder = [ [dA0] , [dA1] , [dA2] ]

END
