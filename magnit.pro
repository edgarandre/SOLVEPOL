pro magnit,image,wcsimage,file=pairfile,cat=catalog
;+
; NAME:
;
;      MAGNIT
;
; PURPOSE:
;
;      PHOTOMETRY AND CALIBRATES THE MEASURED FLUX OF THE STARS 
;      
; EXPLANATION:
;        
;	PERFORMS PHOTOMETRY WITHIN A 2*FWHM DIAMETER APERTURE, AND
;	CALIBRATES THE ZERO POINT MAGNITUD.
;	
;
; CALLING SEQUENCE:
;
;     magnit,image,wcsimage,file=calcpolout,cat=catalog
;
; INPUTS:
;       image - image to perform the photometry
;
;	wcsimage - image with  the wcs information to get ra and dec
;
;       file - file with ID, X, Y and best apertures (from calcpol) organized in columns
;
;       cat - catalog to estimate magzeropoint: VMAG , E_VMAG , RAJ2000 , DEJ2000 
;
;
; OPTIONAL INPUTS:
;     
;
; OUTPUTS:
;     
;
; EXAMPLE:
;
;       magnit,'bias_flat_object0001.fits','fakesky.new.fits',file='bias_flat_object0001.pair',cat='solvepol/GSC2_3.txt'
;
;       
; PROCEDURES USED:
;       
; NOTES:
;
;       The file .pair must exist.
;
; REVISON HISTORY:
;       
;-  	


  subdir='solvepol'
  COMMON FWHMaper,calcfwhm,fwhm

;==values used to run aper (used in runphot)==
  phpadu=5.			
  sky_annulus=[10,20]
  badpix=[-30000,60000]

  img=readfits(image,hdr)

;;===reading stars from .pair file to calculate z_mag===
 filename=subdir+'/'+ pairfile
 readcol,filename,F='I,D,D,D,D,D',starnu,xvalp,yvalp,fluxesp,sharpnp,roundnp,COUNT=ngoodpair,NLINES=nlinespair,/silent

  xvalpo=fltarr(ngoodpair/2)
  yvalpo=fltarr(ngoodpair/2)
  xvalpe=fltarr(ngoodpair/2)
  yvalpe=fltarr(ngoodpair/2)
  magp_oe=fltarr(ngoodpair/2)
  idpo=fltarr(ngoodpair/2)


;apertures in pixels to perform the photometry (radii)
;Good aperture diameter 2 times the fwhm; else r=10.
  if fwhm gt 10 then app_o=10 else app_o=fwhm
  app_e=app_o

  print,'Aperture to perform the photometry, r = ', app_e

;;photometry on ordinay (even) positions
  for i=0, ngoodpair/2-1 do begin
     xvalpo[i]=xvalp[2*i]
     yvalpo[i]=yvalp[2*i]
     idpo[i]=starnu[2*i]
    
     aper,img,xvalpo[i],yvalpo[i],magpo,errmagpo,skypo,skyerropo,phpadu,app_o,sky_annulus,badpix,$
          /SILENT,/flux


;;photometry on extraordinay (pair) positions
     xvalpe[i]=xvalp[2*i+1]
     yvalpe[i]=yvalp[2*i+1]
     ;ide[i]=starnum[2*i+1]

     aper,img,xvalpe[i],yvalpe[i],magpe,errmagpe,skype,skyerrope,phpadu,app_e,sky_annulus,badpix,$
          /SILENT,/flux

;;adding up the magnitude_ordinary +  magnitude_extraordinary
     magp_oe[i] = magpo + magpe
  endfor

;stdev2 = (stdevo2 + stdeve2 ) / 4
;nsky = (nskyo + nskye ) / 2

 ;;calculating the uncalibrated magnitude of all pairs
  zmag_0 = 0.
;  correc = 0.

  EXPTIME=sxpar(hdr,'EXPTIME') ;;read exptime from header
  print,'EXPTIME = ',EXPTIME,' sec'

;==subst comma for point
;  comma = STRPOS(EXPTIME,',')
;  if comma ne -1 then STRPUT, EXPTIME, '.', comma
;  EXPTIME= float(EXPTIME)

  magnitudepair = zmag_0 - 2.5*alog10(magp_oe) ;+ 2.5*alog10(EXPTIME) + correc

;;===calculating postion_ordinary to RA DEC
  ;imagebg=readfits(wcsimage,wcshdr)
  wcshdr=HEADFITS(wcsimage)
  xyad,wcshdr,xvalpo,yvalpo,rap,decp,/CELESTIAL


;;printing temporal file with uncorrected magnitude to be read by zeropoint.pro
  close,1
  openw,1,'tmp_magnit.out'
  printf,1,'ID X Y RA DEC MAG'


  ;;===reading x,y,mag etc from .pair (before P/sigma, chisqrt & fail
  ;;to converge modulation fits)
  print,'Performing photometry to calculate the zero point magnitude'
  print,'using the stars in '+filename

  for i=0, ngoodpair/2-1 do begin
     if FINITE(magnitudepair[i]) eq 1 then begin
        printf,1,format='(I5,2(F8.2," "),F9.5," ",F10.6," ",F9.5)',idpo[i],xvalpo[i],yvalpo[i],rap[i],decp[i],magnitudepair[i]
     endif
  endfor
  close,1

;magnitude_erro = 1.0857 * sqrt((sum-area*msky)/gain+area*stdev2+area2*stdev2/nsky)/(sum-area*msky)

;;function to estimate zmag: zeropoint.pro 
  zmag=zeropoint(magfile='tmp_magnit.out',catfile=catalog)
  spawn,'rm tmp_magnit.out'     ; rm temporal file with uncorrected photometry by zeropoint
  print,'mag zero point =', zmag


;;=recalculating magntiall mags with zmag
  magnitudepair = zmag + magnitudepair ;2.5*alog10(magp_oe) + 2.5*alog10(EXPTIME) + correc

  close,3
  openw,3,subdir+'/'+'magnitall.out'
  printf,3,'zeropointmagnitud=',zmag
  printf,3,'   ID 	X 	Y 	RA 	DEC 	MAG'

  for i=0, ngoodpair/2-1 do begin
     IF Finite(magnitudepair[i]) NE 0  THEN BEGIN

        printf,3,format='(I5,2(F8.2," "),F9.5," ",F10.6," ",F9.5)',idpo[i],xvalpo[i],yvalpo[i],rap[i],decp[i],magnitudepair[i]
       ; print,fix(idpo[i]),rap[i],decp[i],magnitudepair[i]
     endif
  endfor

  print,''
  print,'Calibrated magnitudes of all stars: '+subdir+'/magnitall.out'
  close,3

end
