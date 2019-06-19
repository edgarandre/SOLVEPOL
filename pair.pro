pro pair,filename,shiftx,shifty,deltax,deltay,deltamag
;+
; NAME:
;      PAIR
; PURPOSE:
;      
; EXPLANATION:
;        
;
; CALLING SEQUENCE:
;     PAIR, filename, shiftx, shifty, deltax, deltay, deltamag, [ /MONITOR,
;                       /VINCUT ]
; INPUTS:
;     filename - name of input file containing results of find <filename>.find,
;                or <filename>.cut if vignetting_cut has been run
;     shiftx - shift of two polarization components along the x-axis in pixels
;     shifty - shift of two polarization components along the y-axis in pixels
;     deltax - the +- x-axis search range
;     deltay - the +- y-axis search range
;     deltamag - the magnitude difference limit for finding pairs
;
; OPTIONAL INPUTS:
;     
;
; OPTIONAL KEYWORD INPUTS:
;     
;
; OUTPUTS:
;     
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

subdir = 'solvepol' 
COMMON flux, starnu, fluxes, ngoodpair


ext='.cut'
readcol,filename+ext,F='I,D,D,D,D,D',starnum,xval,yval,flux,sharpness,roundness,COUNT=ngood,NLINES=nlines,/silent

header=strarr(nlines-ngood)

openr,1,filename+ext
readf,1,header
close,1

mag=25-2.5*alog10(flux)
pairnum=0
nopairs=0
nmagcut=0

openw,1,filename+'.pair'
printf,1,transpose(header[0:n_elements(header)-2])
printf,1,''
printf,1,' Program: PAIR '+systime()
printf,1,''
printf,1,' Pairing criteria: shiftx = ',strtrim(shiftx,2),'  shifty = ',strtrim(shifty,2)
printf,1,'     deltax = ',strtrim(deltax,2),'  deltay = ',strtrim(deltay,2),'  deltamag = ',strtrim(deltamag,2)
printf,1,''
printf,1,' No. of sources rejected by SPACIAL    criteria '
printf,1,' No. of sources rejected by DELTAMAG   criteria '
printf,1,''
printf,1,header[n_elements(header)-1]

format='(12x,i5,2f8.2,f9.1,2f9.2)'

for i=0,ngood-1 do begin
	alpha = xval[i]+shiftx
	beta = yval[i]+shifty
	valid_points=where(xval gt alpha-deltax and xval lt alpha+deltax and yval gt beta-deltay and yval lt beta+deltay)
	if valid_points[0] ne -1 then begin
		 valid_mag_points=where(abs(mag(valid_points)-mag(i)) lt deltamag)
		 if valid_mag_points[0] ne -1 then begin
			k=valid_points[valid_mag_points]
			smallest_distance = min((xval[k]-alpha)^2+(yval[k]-beta)^2,best_point)
			pairnum++
			j=k[best_point]

			printf,1,'PAIR',pairnum
			printf,1,f=format,starnum[i],xval[i],yval[i],flux[i],sharpness[i],roundness[i]
			printf,1,f=format,starnum[j],xval[j],yval[j],flux[j],sharpness[j],roundness[j]
		endif else nmagcut++
	endif else nopairs++
endfor
close,1


print,'No. of pairs = ',pairnum
print,'No. of non pairs = ',nopairs

if pairnum gt 32767 then begin
   print,'The number of pairs is greater than the maximum allowed (32767)'
   retall
endif

if pairnum eq 0 then begin
   print,'There are no sources (pairs) left.'
   print,'Verify x- and y-shift, and/or try with lower flux/sigma.'
   print,'...Quitting.'
   retall
endif

n=FILE_LINES(filename+'.pair')
tmp=strarr(1,n)

openr,1,filename+'.pair'
readf,1,tmp
close,1

tmpnew=repstr(tmp,' No. of sources rejected by SPACIAL    criteria ',' No. of sources rejected by SPACIAL    criteria   '+strtrim(nopairs))
tmpnew=repstr(tmpnew,' No. of sources rejected by DELTAMAG   criteria ',' No. of sources rejected by DELTAMAG   criteria   '+strtrim(nmagcut))
openw,1,filename+'.pair'
printf,1,tmpnew
close,1

;;Reading the flux estimated by find and ordered by pair (gaussian flux. This is not the
;;same to the aperture photometry used to calculate the polarisation)
readcol,filename+'.pair',F='I,D,D,D,D,D',starnu,xvalu,yvalu,fluxes,sharpn,roundn,COUNT=ngoodpair,NLINES=nlinespair,/silent ;reading to pass variables

end
