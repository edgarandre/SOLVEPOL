function findsep,filename,x_shift,y_shift
;+
; NAME:
;      FINDSEP
; PURPOSE:
;      
; EXPLANATION:
;        
;
; CALLING SEQUENCE:
;     FINDSEP, filename, x_shift, y_shift
; INPUTS:
;     FILENAME - image array = readfits( <filename>.fits)
;
; OPTIONAL INPUTS:
;    
;
; OPTIONAL KEYWORD INPUTS:
;     
;
; OUTPUTS:
;     X_SHIFT - x-axis shift between ordinary and extraordinary components in pixels
;     Y_SHIFT - y-axis shift between ordinary and extraordinary components in pixels
;
; EXAMPLE:
;       
;       
; PROCEDURES USED:
;       
; NOTES:
;       
; REVISON HISTORY:
;       
;-


  a=filename;readfits(filename,ha)
  asize=size(a)
  border = asize[1]*0.09
  shiftnumber = 7

;Could be asize/4,of asize/9
  result=correl_images(a[border:asize[1]-border-1,border:asize[2]-border-1],a[border:asize[1]-border-1,border:asize[2]-border-1],xshift=asize[1]/4,yshift=asize[2]/4,reduction=8)
  ;window,11
  ;surface,result
  sr=size(result)
  order=sort(-result)               ;max result index
;  ncol=sr(1)                        ;number of columns
;  col=order[0:2] MOD ncol           ;remainder of 3 max / ncolumns; col_index of 3 max
;  row=order[0:2] / ncol             ;3 max / ncolumns; row_index of 3 max
;  x_int_shift=( (col(2)-col(1))/2.) ;half distance in column of matrix between the 2nd and 3rd max 
;  y_int_shift=( (row(2)-row(1))/2.) ;half distance in row of matrix between the 2nd and 3rd max 

  ind = ARRAY_INDICES(result, order[1:2]) 
  x_int_shift=(ind[2]-ind[0])/2.
  y_int_shift=(ind[3]-ind[1])/2.

  print,'Initial estimate of pair separation is :'
  print,'       x shift = ',x_int_shift*shiftnumber
  print,'       y shift = ',y_int_shift*shiftnumber

  result=correl_images(a[border:asize[1]-border-1,border:asize[2]-border-1],a[border:asize[1]-border-1,border:asize[2]-border-1],xoffset_b=x_int_shift*shiftnumber,yoffset_b=y_int_shift*shiftnumber) 
  ;window,11
  ;surface,result
  junk=max(result,maxval)
  sr=size(result)
  ncol=sr(1)
  nrow=sr(2)
  col=maxval MOD ncol                                       ;col_index of max
  row=maxval / ncol                                         ;row_index of max
  x_int_shift=x_int_shift*shiftnumber+(col-((ncol-1)/2))    ;+ distance in column between the max and center 
  y_int_shift=y_int_shift*shiftnumber+(row-((nrow-1)/2))    ;+ distance in row between the max and center


  print,'Closer estimate of pair separation is :'
  print,'       x shift = ',x_int_shift
  print,'       y shift = ',y_int_shift


  result=correl_images(a[border:asize[1]-border-1,border:asize[2]-border-1],a[border:asize[1]-border-1,border:asize[2]-border-1],xoffset_b=x_int_shift,yoffset_b=y_int_shift)
  gauss_result=gauss2dfit(result,gauss_para)
  fwhm=sqrt(2.*alog(2.))*(gauss_para(2)+gauss_para(3))
  ;window,11
  ;surface,result


  cntrd, result, (ncol-1)/2, (nrow-1)/2, xcen, ycen,fwhm; ncol/2
  x_cntrd_shift=x_int_shift+xcen-(ncol-1)/2 ; + distance in column between max and center
  y_cntrd_shift=y_int_shift+ycen-(nrow-1)/2 ; + distance in row between max and center



  if x_cntrd_shift lt 0 and y_cntrd_shift lt 0 then begin
     x_cntrd_shift = sqrt(x_cntrd_shift*x_cntrd_shift)
     y_cntrd_shift = sqrt(y_cntrd_shift*y_cntrd_shift)
  endif

  if x_cntrd_shift lt 0 and y_cntrd_shift gt 0 then begin
     x_cntrd_shift = sqrt(x_cntrd_shift*x_cntrd_shift)
     y_cntrd_shift = - sqrt(y_cntrd_shift*y_cntrd_shift)
  endif

  
  print,'Final pair separation after centroiding is:'
  print,'       x shift = ',x_cntrd_shift
  print,'       y shift = ',y_cntrd_shift


     
  x_shift=x_cntrd_shift
  y_shift=y_cntrd_shift

  
  return, result

end
