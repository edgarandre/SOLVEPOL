pro astrometrylist,image,calcpol

;+
; NAME:
;      ASTROMETRYLIST
; PURPOSE:
;      This procedure puts the wcs to the given image (fakesky.fits). 
;           
;      
; EXPLANATION: 
;       It performes aperture photometry of the ordinary stars given
;       in calcpol. A text file with the position (Xo,Yo) of stars
;       sorted by the brigthes is printed. Enter the text file to the
;       Astrometry.net site and download the wcs.fits from the results
;       page, and give it to this procedure. The output is
;       fakesky.new.list.fits with the wcs information.
;
;
;
; CALLING SEQUENCE:
;
;     ASTROMETRYLIST,image,calcpol
;
; INPUTS:
;
;       image - fits image to proccess. Preferently fakesky.fits
;
;       calcpol - file text output of calcpol listing the 
;
; OPTIONAL INPUTS:
;     
;      
;
;
;
; OUTPUTS:
;
;      textXYsorted.list - text file listing the position X and Y in
;                          pixels of the brigthest stars in decressing
;                          order. 
;
;      fakesky.new.list.fits - Image with thte wcs information.
;
; EXAMPLE:
;       
;       ASTROMETRYLIST,'fakesky.fits','~/Desktop/hd110v/solvepol/calcpol.out'
;
; PROCEDURES USED:
;       
; NOTES:
;       
; REVISON HISTORY:
;       
;-  	
;	

;;===reading image===
  img=readfits(image,hdr)


;;===Detected stars in calcpol===
  calcpolout=calcpol
  readcol,calcpolout,format='I,D,D,I',starnum,xval,yval,bestaprid,COUNT=ngood,NLINES=nlines,/silent


  phpadu=5			;values used to run aper (used in runphot of calcpol)
  sky_annulus=[10,20]
  badpix=[-30000,60000]

;===position X and Y image=== 
;reading Xo and Yo ordinary (even)
  Xo=fltarr(ngood/2)
  Yo=fltarr(ngood/2)
  ido=fltarr(ngood/2)

  for i=0, ngood/2-1 do begin
     Xo[i]=xval[2*i]
     Yo[i]=yval[2*i]
     ido[i]=starnum[2*i]
  endfor

;;photometry on ordinay (even) positions
  apper=7
  aper,img,Xo,Yo,magso,errmagso,skyo,skyerro,phpadu,apper,sky_annulus,badpix,/flux,/silent

;sort (magso) by counts
  magsosortID=bsort(magso,magsosort,/info,/reverse)

;;open a file listing the X,Y positions of sources, sorted with
;;the brightest sources firs
  output='textXYsorted.list'
  close,3
  openw,3,output
  for i=0,ngood/2-1 do begin
     printf,3,F='(F8.2,F8.2)',Xo[magsosortID[i]],Yo[magsosortID[i]]
  endfor
  close,3


  print,'astrometrylist output: textXYsorted.list'
;  print,''
;  wcsfile=''
;  read,'Give me the wcs.fits (output of Astrometry.net):',wcsfile

;;to put the wcs into the header
;  img2=readfits(wcsfile,hdrwcs)           ;output of Astrometry.net
;  extast,hdrwcs,astr                      ;extract the astrometric information 
;  putast,hdr,astr                         ;put the astrometry into the hdr
;  writefits,'fakesky.new.list.fits',img,hdr ;write the image with the complete wcs info.

;  print,'===Output fakesky.new.list.fits==='


end
