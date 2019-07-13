pro findstdstar,pairfile,wcsimage,stdstars,RA,DEC
;+
; NAME:
;
;      FINDSTDSTAR
;
; PURPOSE:
;
;      FINDS STANDARD STAR 
;      
; EXPLANATION:
;        
;	WHEN THE FIELD HAS A STANDARD STAR, THE MEASURED THETA AND
;	THETA FROM THE LITERATURE ARE COMPARED TO GET DELTATHETA
;	CORRECTION.
;	
;
; CALLING SEQUENCE:
;
;     deltatheta,fintab,polstdstars
;
; INPUTS:
;
;       pairfile - .pair text file.
;
;       wcsimage - fits file with the wcs information in the header
;
;       stdstars - FILE WITH INFORMATION FROM THE LITERATURE OF THE
;                     STANDARDSTARS: polarimetrystandardstars.dat
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
;      FINDSTDSTAR,'solvepol/bias_flat_hd155197v0001.pair','wcs.fits','/Users/Edgar/Desktop/pipeline/polstdstars.dat',RA,DEC 
;
;       
; PROCEDURES USED:
;       
; NOTES:
;       
; REVISON HISTORY:
;       
;-  	

;;===reading stars from .pair file to find standard star===
  readcol,pairfile,F='I,F,F,D,D,D',starnu,xvalp,yvalp,fluxesp,sharpnp,roundnp,COUNT=ngoodpair,NLINES=nlinespair,/silent
  readcol,pairfile,F='A,I',Pair,Pairno,/silent

  xvalpo=fltarr(ngoodpair/2)
  yvalpo=fltarr(ngoodpair/2)
  xvalpe=fltarr(ngoodpair/2)
  yvalpe=fltarr(ngoodpair/2)
  idpo=intarr(ngoodpair/2)
  idpe=intarr(ngoodpair/2)


  for i=0, ngoodpair/2-1 do begin
     ;;ordinay (even) positions
     xvalpo[i]=xvalp[2*i]
     yvalpo[i]=yvalp[2*i]
     idpo[i]=starnu[2*i]

     ;;extraordinay (pair) positions
     xvalpe[i]=xvalp[2*i+1]
     yvalpe[i]=yvalp[2*i+1]
     idpe[i]=starnu[2*i+1]
  endfor


;===calculating postion_ordinary to RA DEC
  ;imagebg=readfits(wcsimage,wcshdr)
  wcshdr=HEADFITS(wcsimage)
  xyad,wcshdr,xvalpo,yvalpo,rap,decp,/CELESTIAL


;==Read literature catalogue.
  readcol,stdstars,format='A,D,D,D,D,D,A',name,RAstd,DECstd,thetastd,pstd,dpstd,Refstd,COUNT=ngoodstd,NLINES=nlinesstd,/silent


;==match RA & DEC
  sigmaRA=0.001                ;~3.6 degrees 
  sigmaDEC=0.001 

  starfind=0

  for i = 0,ngoodpair/2-1 do begin
     for j=0,ngoodstd-1 do begin

        if (Rap[i] ge (Rastd[j]-sigmaRA)) and (Rap[i] le (Rastd[j]+sigmaRA)) $
           and (Decp[i] ge (Decstd[j]-sigmaDEC)) and (Decp[i] le (Decstd[j]+sigmaDEC)) then begin

           ;==if match then take theta & theta_cat
           print,' --- ---'
           print,'Standard star on the field: ',name[j]
           print,'RA & DEC:',Rap[i],Decp[i]

           print,pair[i*3],pairno[i*3]
           pair= pair[i*3]
           pairno = pairno[i*3]

           print,'   STAR       X       Y'
           print,idpo[i], xvalpo[i],yvalpo[i],FORMAT='(I5," ",2(F7.2," "))'
           ido=idpo[i]
           xo=xvalpo[i]
           yo=yvalpo[i]
           print,idpe[i], xvalpe[i],yvalpe[i],FORMAT='(I5," ",2(F7.2," "))'

           print,''

           pstdl=pstd[j]
           dpstdl=dpstd[j]
           thetastdl=thetastd[j]
           Refstdl=Refstd[j]

           starfind=1

        endif 

     endfor
  endfor

  if starfind eq 0 then begin
     print,'There is not a standard star on the FOV.' 
     return
  endif


;===reading <name>all.pol === 
  lengthname=strlen(pairfile)
  namesdim=lengthname-5
  polfile= strmid(pairfile,0,namesdim)+'all.pol'


  close,1
  openr,1,polfile
  nlines=file_lines(polfile)
  tmpstr=strarr(nlines)
  readf,1,tmpstr
  close,1

  ;print,'nlines:',nlines
  for i=0L,nlines-2 do begin
                                ;  ;uncomente/change if there is an
                                ;  error here.
  ;for i=0,nlines-2 do begin
     if strcmp(tmpstr[i],Pair+STRTRIM(pairno)) then k=i
  endfor

;==finding index of apperture
  indexapper=strpos(tmpstr[k+64],'psigma):')
  tmpapper=float(strmid(tmpstr[k+64],indexapper+9))


;==reading measured z & dz 
  print,'Best Aperture (aperture with lowest psigma):',fix(tmpapper)

  print,tmpstr[k+3  +6*fix(tmpapper)-6+3]
  print,tmpstr[k+3  +6*fix(tmpapper)-6+4]
      
  zs=STRSPLIT(tmpstr[k+3  +6*fix(tmpapper)-6+3],/EXTRACT)
  dzs=STRSPLIT(tmpstr[k+3  +6*fix(tmpapper)-6+4],/EXTRACT)

  z=float(zs[2:*])
  dz=float(dzs[2:*])

;===printing modulation
  print,tmpstr[k+3  +6*fix(tmpapper)-6+5]
  print,tmpstr[k+3  +6*fix(tmpapper)-6+6]
  print,''


;==calculating Chisqrt
  chisqrt=graf(z,dz,/plotm)
  print,'Reduced Chi-square=',chisqrt

;===printing polarimetric properties
  strpolprop=STRSPLIT(tmpstr[k+3  +6*fix(tmpapper)-6+6],/EXTRACT)
  p=float(strpolprop[4])
  dp=float(strpolprop[5])
  theta=float(strpolprop[7])
  dtheta=float(strpolprop[8])

  print,'' 
  print,'Pol_measured & Theta_measured =',p*100.,'+/-',dp*100.,'%, ', theta,'+/-',dtheta,' degrees', FORMAT='(4(A,F7.3),A)'
  print,'Pol_literature & Theta_literature = ',pstdl,'+/-',dpstdl,'%, ',thetastdl,' degrees '+Refstdl,FORMAT='(A,F7.3,A,F7.3,A,F7.3,A)' 

  ;==print theta_literature - theta_measured
  deltatheta =  thetastdl -  theta 
  print,'Deltatheta (degrees) = theta_literature - theta_measured = ',deltatheta 


  fintabtrue=0
  
  ngood=n_elements(RA)
  for i = 0,ngood-1 do begin
     for j=0,ngoodstd-1 do begin

        if (Ra[i] ge (Rastd[j]-sigmaRA)) and (Ra[i] le (Rastd[j]+sigmaRA)) $
           and (Dec[i] ge (Decstd[j]-sigmaDEC)) and (Dec[i] le (Decstd[j]+sigmaDEC)) then begin

           fintabtrue=1
           print,'Standard star on the final catalogue: ',name[j]

        endif 
     endfor
  endfor

  if starfind eq 0 and fintabtrue eq 0 then print,'There is not a standard star on the FOV.'

  if starfind eq 1 and fintabtrue eq 0 then print,'The standard star is not in the final catalogue.'

  if starfind eq 1 and fintabtrue eq 1 then print,'The standard star is on the FOV and in the final catalogue.'


  print,' --- ---'

end
