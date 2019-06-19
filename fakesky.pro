pro fakesky,image,pairfile
;+
; NAME:
;      FAKESKY
; PURPOSE:
;      MASKOUT THE STARS ON THE EXTRAORDINARY STARS.
;      THE PROCEDURE CAN BE MODIFY TO MASKOUT THE ORDINARY STARS.
;      PAIRFILE WITH POSITIONS IS OUTPUT OF POLPAIR.PRO.
;      
; EXPLANATION:
;        setup path 
;	;!PATH = Expand_Path('+/Users/Edgar/IDLWorkspace80/Default/pipeline/') + ':' + !PATH
;	
;
; CALLING SEQUENCE:
;     FAKESKYV3, image, pairfile
; INPUTS:
;       image - image of reference to produce  the fake sky
;
;	pairfile - file with positions (output of polpair.pro) in
;                  columns: ID   X       Y       BestApp	 
;
; OPTIONAL INPUTS:
;     
;
; OUTPUTS:
;     
;
; EXAMPLE:
;
;       FAKESKY,'hd110cv0001.fits','calcpol.out'
;       
; PROCEDURES USED:
;       
; NOTES:
;       
; REVISON HISTORY:
;       
;- 


;read image===
  image1=readfits(image,header) ;Reading fits file


;===position X and Y  image===  
  readcol,pairfile,F='I,D,D,I',starnum,xval,yval,BestApp,COUNT=ngood,NLINES=nlines,/silent


  if ngood eq 0 then begin
     print,'There are no sources (pairs) left.'
     print,'Verify x- and y-shift, or try with lower flux/sigma and/or pol/sigma_pol.'
     print,'...Quitting.'
     retall
  endif


;redifine X and Y
  X=fltarr(ngood/2)
  Y=fltarr(ngood/2)
  radius=fltarr(ngood/2)

;for i=0, ngood/2-1 do begin ;ORDINARY
;X[i]=xval[2*i]
;Y[i]=yval[2*i]
;endfor

  for i=0, ngood/2-1 do begin   ;EXTRAORDINARY
     X[i]=xval[2*i+1]
     Y[i]=yval[2*i+1]
     radius[i]=BestApp[2*i+1]
  endfor


;===mode of the sky====
  sky,image1,skymode,skysig
  value=skymode

;======Making fake background====
  sizex = n_elements(image1[*,0])
  sizey = n_elements(image1[0,*])
  map_size = size(image1)


  newmap = fltarr(map_size[1]+3*10,map_size[2]+3*10)
  newmap[15:map_size[1]+15-1,15:map_size[2]+15-1] =  image1

  for j=0,ngood/2-1 do begin    ;EXTRAORDINARY
     
     xint = round(X[j])
     yint = round(Y[j])
     rint = fix(radius[j])
  
     maskarray=fltarr(rint*3,rint*3)
     maskarray[*,*]=value
     maskarray[*,*] = maskarray[*,*] + RANDOMN(skymode,rint*3,rint*3)

     maskarray2=fltarr(rint,rint)
     maskarray2[*,*]=value
     maskarray2[*,*] = maskarray2[*,*] + RANDOMN(skymode,rint,rint)
     

 
;==Masking===
     ;;square mask with 3r side
     newmap[xint+15-1.5*rint:xint+15+1.5*rint-1 , yint+15-1.5*rint:yint+15+1.5*rint-1] = maskarray 

     ;;square mask with r to make a circular mask 
     newmap[xint+15-rint*2:xint+15-1-rint*1 , yint+15-rint*0.5:yint+15-1+rint*0.5] = maskarray2
     newmap[xint+15+rint*1:xint+15-1+rint*2 , yint+15-rint*0.5:yint+15-1+rint*0.5] = maskarray2

     newmap[xint+15-rint*0.5:xint+15-1+rint*0.5 , yint+15-rint*2:yint+15-1-rint*1] = maskarray2
     newmap[xint+15-rint*0.5:xint+15-1+rint*0.5 , yint+15+rint*1:yint+15-1+rint*2] = maskarray2 


  endfor


  newmap = newmap[15:map_size[1]+15-1,15:map_size[2]+15-1]
  writefits,'fakesky.fits',newmap,header


  print,''
  print,'FAKESKY run well: output fakesky.fits'
  print,''


end
