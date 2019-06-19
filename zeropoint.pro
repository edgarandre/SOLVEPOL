function zeropoint,magfile=magfile,catfile=catfile
;+
; NAME:
;
;      ZEROPOINT
;
; PURPOSE:
;
;      CALCULATES THE ZERO POINT MAGNITUDE 
;      
; EXPLANATION:
;        
;	RETURNS THE ZEROPOINT MAGNITUDE
;	
;
; CALLING SEQUENCE:
;
;     zeropoint,magfile=magfile,catfile=catfile
;
; INPUTS:
;       magfile - File (temporal: tmp_magnit.out) with
;       ido[i],Xo[i],Yo[i],ra[i],dec[i],magnitude[i]. Where Magnitud
;       is not corrected by zero point
;
;
;       catfile - catalog with Format=D,D,D,D (Mag,sigmaMAg,Ra,Dec)
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
;       zmag=zeropoint(magfile='tmp_magnit.out',catfile=catalog)
;
;       
; PROCEDURES USED:
;       
; NOTES:
;       
; REVISON HISTORY:
;       
;-  	


  readcol,magfile,F='I,D,D,D,D,D',ID,Xo,Yo,Ra,Dec,mydata,COUNT=ngood,NLINES=nlines,/silent 

;magfile='tmp_magnit.out'
;close,1
;openr,1,magfile
;nolines= FILE_LINES(magfile)
;mydata=fltarr(6,nolines) ;this version has problems to read magfile.out with lines *******.
;readf,1,mydata
;close,1

  readcol,catfile,F='D,D,D,D',catMag,catMagerr,catRa,catDec,COUNT=ngoodCat,NLINES=nlinesCat,/silent,COMMENT='N'

  sigmaRA=0.001
  sigmaDEC=0.001

  close,3
  openw,3,'tmp_mag.txt'
  for i = 0,ngood-1 do begin
     for j=0,ngoodCat-1 do begin
;if (mydata[3,i] ge (cat[1,j]-sigmaRA)) and (mydata[3,i] le (cat[1,j]+sigmaRA)) $
;   and (mydata[4,i] ge (cat[2,j]-sigmaDEC)) and (mydata[4,i] le (cat[2,j]+sigmaDEC)) then begin
;mags=[mydata[5,i],cat[3,j]]

        if (Ra[i] ge (catRa[j]-sigmaRA)) and (Ra[i] le (catRa[j]+sigmaRA)) $
           and (Dec[i] ge (catDec[j]-sigmaDEC)) and (Dec[i] le (catDec[j]+sigmaDEC)) then begin
   
           mags=[ID[i],mydata[i],catMag[j]]
           printf,3,mags,f='(I5," ",2(f9.5," "))'

        endif
     endfor
  endfor
  close,3

;===checking if there are matching stars==
  readcol,'tmp_mag.txt',F='I,D,D',ID,mymag,catmag,COUNT=ngoodtmp,NLINES=nlinestmp,/silent 

  if (nlinestmp eq 0) then begin
     print,''
     print,"There are no matching stars between the catalogue and photometry to calibrate the magnitud."

     zero = 0.0
     print,'Using mag zero point = ', zero
     print, ''

     spawn,'rm tmp_mag.txt'     ;temporal file containing matched magnitudes

     return,zero
  endif 

;;sort magnitudes in increasing order (lower mag = brighter star)
  mymagsortID=bsort(mymag,mymagsort,/info)

;;mean of the difference of magnitude of the first 25 brighter stars (Nbrightest_stars)

  if (ngoodtmp le 25) then begin
     Nbrightest_stars=ngoodtmp
  endif else begin
     Nbrightest_stars=25
  endelse

  print,'Brightest stars (Nbrightest_stars) = ',Nbrightest_stars
  print,'ID:',ID[0:Nbrightest_stars-1]

  DeltaMag=0

  for j=0,Nbrightest_stars-1 do begin
     DeltaMag=DeltaMag+ABS((catmag[mymagsortID[j]]-mymagsort[j]))
     zero=DeltaMag/Nbrightest_stars
  endfor

  spawn,'rm tmp_mag.txt'        ;temporal file containing matched magnitudes

;  print,'mag zero point =',zero
  print,'==== zeropoint succesfully finished.  ==='
  return,zero
end
