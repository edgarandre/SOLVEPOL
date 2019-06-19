pro polpair,objectfile,NWAVEPLATE=nwaveplate,deltatheta=deltatheta
;+
; NAME:
;      POLPAIR
; PURPOSE:
;      
; EXPLANATION:
;        
;
; CALLING SEQUENCE:
;     POLPAIR, objectfile, nwaveplate=
; INPUTS:
;     objectfile - text file listing image names, without file extensions
;                  (i.e., without '.fits'), order of waveplate position
;     nwaveplate - number of waveplate positions
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

  subdir='solvepol'
  COMMON sig,skysig
  COMMON FWHMaper,calcfwhm,fwhm
  COMMON Fluxsigma,numsig

;if N_params() eq 0 then begin
;      print,'Syntax - ' + $
;        'polpair, objectfile, NWAVEPLATE=, VIGNCUT=, deltatheta=' 
;      return
;endif


;if not keyword_set(NWAVEPLATE) then begin
;    nwaveplate=''
;    read,nwaveplate,prompt='How many waveplate positions? (4, 8, or type Q to quit)'
;    while nwaveplate ne '4' and nwaveplate ne '8' and nwaveplate ne 'Q' and n;waveplate ne 'q' do begin
;	print,'The value entered is not an accepted value. Please try again.'
;	read,nwaveplate,prompt='How many waveplate positions? (4, 8, or type Q to quit)'
;    endwhile
;    if nwaveplate eq 'Q' or nwaveplate eq 'q' then return
;endif

; if not keyword_set(ROTATED) then rotated='no'

  nobjfiles = FILE_LINES(objectfile) ;List of files group by field in the order of waveplate postion
  objfiles=strarr(nobjfiles)

  openr,1,objectfile
  readf,1,objfiles
  close,1


;removing the .fits extenction to get the name only
  names=strarr(nobjfiles)
  for i=0,nobjfiles-1 do names[i]=strmid(objfiles[i],0,strlen(objfiles[i])-5)


;***********paramters that can be adapted/changed********************************

;numsig: the number of standard deviation above the background for detection. Used by daofind
;numsig = ''
;read,'Give me flux/sigma value : ', numsig		

  sharplim = [-5.0,5.0]         ;the sharpness limits for the find routine (the default listed for find is [0.2,1.0])
  roundlim = [-1.5,1.5]         ;the roundness limits for the find routine (the default listed for find is [-1.0,1.0])

  deltax=2                      ;values used by pair routine
  deltay=2
  deltamag=10

  ;-values used by runphot routine (passed to aper)-
  phpadu=5.			;Photons per Analog per Digital Unit
  apr=[1,2,3,4,5,6,7,8,9,10]
  skyrad=[10,20]
  badpix=[-30000,60000]

;rdn - value used by calcpol routine
  img=readfits(objfiles[0],hdr)
  rdn=sxpar(hdr,'RDNOISE') ;;read rdn from header
  print,'RDNOISE = ',rdn

;*******************************************

  cur_wpp=1                     ;This is the current waveplate position.

  for i=0,nobjfiles-1 do begin  ;start a loop to run daofind, ordem and phot on all images

     tmpimg=readfits(objfiles[i],tmph)

     sky,tmpimg,skymode,skysig

     print,'For image: ',objfiles[i],' the sky background has a mode of ',skymode,' and a sigma of ',skysig

 
;determine the point source locations with find, and pair them.
;These coords will be used as reference for all other positions of the instrument.
     if cur_wpp eq 1 then begin	;First waveplate position only			


        refimg=readfits(objfiles[i],hrefimg)
        size=size(refimg)

        if calcfwhm eq 1 then begin
           ;fwhm = GAUSSFIT_INPARTS(objfiles[i])
           calcufwhm,objfiles[i],fwhm
           print,''
           print,'The average FWHM of the reference image is approximately: ',strtrim(fwhm,2)
        endif else print,'The FWHM set to: ',strtrim(fwhm,2)



        ;numsig: the number of standard deviation above the background for detection. Used by daofind
        print,'flux/sigma value : ',numsig

	;run find to find stars in image and write to <object>.find
        hmin = numsig*skysig
        find, tmpimg, x, y, flux, sharp, roundness, hmin, fwhm, roundlim, sharplim, print=subdir+'/'+names[i]+'.find'


        ;run findsep to find separation between ord. and ext.==
        print,'Run findsep to find the separation of the ordinary and extraordinary components'

        find, tmpimg, xsep, ysep, fluxsep, sharpsep, roundnesssep, 10*skysig, fwhm, roundlim, sharplim ;to findsep
        sizex = n_elements(tmpimg[*,0])
        sizey = n_elements(tmpimg[0,*])
        tmpimage=fltarr(sizex+6*2,sizey+6*2)
        maskarray=fltarr(6,6) 
        maskarray[*,*]=5

        for j=0,n_elements(xsep)-1 do begin
           tmpimage[6+xsep[j]:6+xsep[j]+6-1 , 6+ysep[j]:6+ysep[j]+6-1] = maskarray ;==to findsep
        endfor

;        writefits,'tmpimage.fits',tmpimage ;==to findsep        
        correl_result=findsep(tmpimage,shiftx,shifty)

        if shiftx eq 0. and shifty eq 0. or ABS(shiftx) gt 200 or ABS(shifty) gt 200 then begin
           print,'Second try to find the separation using the reference image'
           filename=readfits(objfiles[i])
           correl_result=findsep(filename,shiftx,shifty)
        endif

        if shiftx eq 0. and shifty eq 0. or ABS(shiftx) gt 200 or ABS(shifty) gt 200 then begin

           answer = ''
           read,answer,prompt='Are x shift and y shift correct? (y/n):'

           if (answer eq 'y') or (answer eq 'yes') then begin
           
              shiftx=shiftx
              shifty=shifty

           endif else begin
              x_cntrd_shift=''
              y_cntrd_shift=''
              read,x_cntrd_shift,prompt='Input x shift:'
              read,y_cntrd_shift,prompt='Input y shift:'
           
              shiftx=x_cntrd_shift
              shifty=y_cntrd_shift
           endelse
        endif
       


        ;run the vignetting on current image to remove detections in the vignetted regions
        vignetting,objfiles[i]			


        ;run pair to find ordinary/extraordinary pairs
        print,'running pair'
        pair,subdir+'/'+names[i],shiftx,shifty,deltax,deltay,deltamag
			
     endif else begin
	 
	;perform cross-correlation to determine image shift from first waveplate position
        imgb=readfits(objfiles[i],himgb)
        sizerefimg=size(refimg)
        sizeimgb=size(imgb)
;        correl_optimize,refimg[20:sizerefimg[1]-20,20:sizerefimg[2]-20],imgb[20:sizeimgb[1]-20,20:sizeimgb[2]-20],xoffset_optimum,yoffset_optimum,/NUMPIX,magnification=2
;        correl_optimize,refimg,imgb,xoffset_optimum,yoffset_optimum,/NUMPIX,magnification=2
        correl_optimize,refimg[sizerefimg[1]*.1:sizerefimg[1]-sizerefimg[1]*.1,sizerefimg[1]*.1:sizerefimg[2]-sizerefimg[1]*.1],imgb[sizeimgb[1]*.1:sizeimgb[1]-sizeimgb[1]*.1,sizeimgb[1]*.1:sizeimgb[2]-sizeimgb[1]*.1],xoffset_optimum,yoffset_optimum,/NUMPIX,magnification=2

        xoffset_optimum = xoffset_optimum
        yoffset_optimum = yoffset_optimum

        print,'For waveplate position ['+STRTRIM(i+1,2)+'] , the offset is:',xoffset_optimum,yoffset_optimum
        readcol,subdir+'/'+names[i-cur_wpp+1]+'.pair',F='L,D,D,D,D,D',starnum,xval,yval,flux,sharpness,roundness,COUNT=ngood,NLINES=nlines,/silent
        header=strarr(nlines-ngood-ngood/2)
        openr,1,subdir+'/'+names[i-cur_wpp+1]+'.pair'
        readf,1,header
        close,1
        format='(12x,i5,2f8.2,f9.1,2f9.2)'
        openw,1,subdir+'/'+names[i]+'.pair'
        printf,1,transpose(header)
        for k=0,ngood/2-1 do begin
           printf,1,'PAIR',k+1
           printf,1,f=format,starnum[k*2],xval[k*2]-xoffset_optimum,yval[k*2]-yoffset_optimum,flux[k*2],sharpness[k*2],roundness[k*2]
           printf,1,f=format,starnum[k*2+1],xval[k*2+1]-xoffset_optimum,yval[k*2+1]-yoffset_optimum,flux[k*2+1],sharpness[k*2+1],roundness[k*2+1]
        endfor
        close,1
     endelse
     print,'run runphot to determine the photometry'
     runphot,names[i],phpadu,apr,skyrad,badpix,fwhm
	
	;on the last waveplate position run calcpol to calculate the polarization
     if cur_wpp eq nwaveplate then begin
        cur_wpp=1               ;reset cur_wpp to 1 for the next field

        
        print, 'running calcpol'
        calcpol,[names[i-nwaveplate+1:i]],nwaveplate,rdn,deltatheta
                        
     endif else cur_wpp++	;increase cur_wpp by 1 for the next position

  endfor

end
