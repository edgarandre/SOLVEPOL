pro SOLVEPOL,objectfile,nreads=nreads,flats=flatfile,bias=biasfile,deltatheta=deltatheta,$
             fluxsigma=fluxsigma,polsigma=polsigma,fwhm=fwhmax,astrometry=astrometry,oplotid=oplotid, $
             Chisqrtt=Chisqrtt,chisqrtint=chisqrtint
;+
; NAME:
;
;      SOLVEPOL
;
; PURPOSE:
;
;      CALCULATES THE POLARIMETRIC, ASTROMETRIC AND PHOTOMETRIC
;      PROPERTIES OF THE DETECTED SOURCES IN THE FIELD
;      SOLVEPOL IS THE MAIN PROCEDURE.
;      
;      
; EXPLANATION:
;
;       HELP in the README.pdf.
;
;       Setup the path to the pipeline folder: 
;
;	!PATH = Expand_Path('+/Users/eramirez/pipelineV16.1/') + ':' + !PATH
;
;       Set in the first lines of these comments the path to the
;       location directory of the polarimetrystdstars.dat file:
;
;       solvepolDir = '+/Users/Edgar/Desktop/pipelineV16.1/'
;
; CALLING SEQUENCE:
;
;     solvepol, objectfile[, nreads = nreads, flats = flatfile, bias =
;     biasfile, deltatheta = deltatheta, fluxsigma = fluxsigma,
;     polsigma = polsigma, fwhm = fwhmax, /astrometry, /oplotid,
;     Chisqrt = Chisqrt, /chisqrtint] 
;
; INPUTS:
;
;       objectfile - text file with the list of images to process
;
; OPTIONAL INPUTS:
;     
;	nreads - the numbebr of readouts per waveplate position
;      
;       flats - text file with the list of flats
;
;       bias - text file with the list of bias
;
;       deltatheta - theta correction
;
;       fluxsigma - flux/sigma_flux value to detect sources. If blank default is 2.
;
;       polsigma - pol/sigma_pol value. If blank default is 1.
;
;       fwhm - The fwhm to find the sources. If blank the pipelien will
;             calculate the average fwhm of the image in the
;             first waveplate position angle. Set the fwhm if you
;             estimate that the fwhm calculated by the pipeline is wrong.
;
;       astrometry - Flag to indicate that the computer has
;                    Astronmetry.net installed. If blank indicates
;                    that the computer does not have Astrometry.net
;                    installed and the pipeline will prompt you to
;                    input the wcs solution manually. Acepted values
;                    for astrometry are 'y' (within quotationmarks) or
;                    /astrometry.
;
;
;       oplotid - flag to over plot the IDs of the stars on the
;                polarimetry vector image. Acepted values are 'y'
;                (within quotationmarks) or /oplotid.
;
;       chisqrt - Limit (number) of the Chisqrt of the fit to the
;                 modulation. Sources with calculated Chisqrt greater
;                 than the given values are rejected. Default Chisqrt
;                 limit is 6.0.
;
;       chisqrtint - Flag to choose the Chisqrt limit in interactive
;                    mode. If blank the pipeline will not stop and the
;                    Chisqrt limit (given or default) is used. Acepted
;                    values are /chisqrtint or 'y' (within
;                    quotationmarks). 
;
; OUTPUTS:
;
;      Explained in README.pdf      
;
; EXAMPLE:
;       
;      IDL> solvepol,'objects.list',deltatheta = 45.0
;       
; PROCEDURES USED:
;       
; NOTES: 
;
;       Read the README.pdf for help
;
; REVISON HISTORY:
;       
;-  	
;	
;Modify the next line to the path to find polarisedstdstars.dat
  solvepolDir = '/Users/eramirez/pipelineV16.1/' 

  COMMON imgnames,name
  COMMON Fluxsigma,numsig
  COMMON Polsigma,p_over_sigma,Chisqrt_lim,chisqrt_int
  COMMON FWHMaper,calcfwhm,fwhm ;polpair,magnit

  subdir='solvepol'
  FILE_MKDIR, subdir
  journal,'solvepol.log'       ;Printing terminal to a log file.
  set_plot,'x'
  !P.MULTI = 0 
 
  if KEYWORD_SET(deltatheta) eq 0 then deltatheta=0.

  if KEYWORD_SET(fluxsigma) eq 0 then numsig=2. else numsig=fluxsigma

  if KEYWORD_SET(polsigma) eq 0 then p_over_sigma=1. else p_over_sigma=polsigma

  if KEYWORD_SET(fwhmax) eq 0 then calcfwhm=1 else begin
     calcfwhm=0 
     fwhm=fwhmax
  endelse

  if KEYWORD_SET(astrometry) eq 1 then astrome='y' else astrome='n'

  if KEYWORD_SET(oplotid) eq 1 then overplotID='y' else overplotID='n'

  if KEYWORD_SET(Chisqrtt) eq 0 then Chisqrt_lim = 6. else Chisqrt_lim = Chisqrtt

  if KEYWORD_SET(chisqrtint) eq 1 then chisqrt_int=1 else chisqrt_int = 0

;==Check if the images have been reduced========
  noImages = FILE_LINES(objectfile) ;Number of images 
  objnames=strarr(noImages)         ;define array with number of images

  close,1
  openr,1,objectfile            ;Open text file with object names
  readf,1,objnames
  close,1

;==Name of file [0] to add as prefix===
  name=strmid(objnames[0],0,strlen(objnames[0])-5)
;========

  tmpimage=readfits(objnames[0],tmphdr)
  tmpimagesize=size(tmpimage)
  print,format='(%"Size of image/array : %i by %i")', tmpimagesize[1],tmpimagesize[2] 

  CCDPROC=sxpar(tmphdr,'CCDPROC')        ;;read  from header
  reducedIDL=sxpar(tmphdr,'Solvepol') ; 'Images corrected by Flats & Bias'


  if (Keyword_Set(CCDPROC) eq 1) xor (Keyword_Set(reducedIDL) eq 1) then begin ;'Images proccessed by reductionv2.pro
     print,'Images already reduced...proceeding with calpol'

     noImages = FILE_LINES(objectfile) ;Number of images 
     objnames=strarr(noImages)         ;define array with number of images

     close,1
     openr,1,objectfile         ;Open text file with object names
     readf,1,objnames
     close,1

     print,'Creating text file with reduced images: objectsreduced.list'
     filewithreducedimages='objectsreduced.list'

     close,2
     openw,2,filewithreducedimages
     for i=0,noImages-1 do printf,2,objnames[i]; strmid(objnames[i],0,namesdim[i])
     close,2

     print,'====processing reduced data===='


  endif else begin

     if (KEYWORD_SET(flatfile) eq 1) or (KEYWORD_SET(biasfile) eq 1) then begin 
        print,'Starting reduction'

        reduction,objectfile,flats=flatfile,bias=biasfile  ;reducing data
        filewithreducedimages=subdir + '/' + 'objectsreduced.list' ;output of reduction

        print,'====data reduced===='

     endif else begin 

        print,'Bias and/or flats files missing: Assuming your data are reduced'

        noImages = FILE_LINES(objectfile) ;Number of images 
        objnames=strarr(noImages)         ;define array with number of images

        close,1
        openr,1,objectfile      ;Open text file with object names
        readf,1,objnames
        close,1

        print,'Creating text file with reduced images: objectsreduced.list'
        filewithreducedimages=subdir + '/' +'objectsreduced.list'

        close,2
        openw,2,filewithreducedimages
        noImages=fix(noimages[0])
        for i=0,noImages-1 do printf,2,objnames[i] 
        close,2
        
        print,'====processing reduced data===='

     endelse                    ;=======

  endelse


;;===combining images with multireads per waveplete position===
  if KEYWORD_SET(nreads) eq 1 then begin 
     dither,filewithreducedimages,dither=nreads

     noImages = FILE_LINES('objectscombined.list' ) ;Number of images 
     objnames=strarr(noImages)                      ;define array with number of images

     spawn,'mv  objectscombined.list '+subdir+'/'
     filewithreducedimages=subdir + '/' +'objectscombined.list'

     print,'Images combined'

  endif


  nReducedFiles=FILE_LINES(filewithreducedimages) ;text file with names of reduced images without .fits
  polpair,filewithreducedimages,NWAVEPLATE=nReducedFiles,deltatheta=deltatheta
  calcpolout=subdir+'/'+'calcpol.out'

  ReducedFiles=strarr(nReducedFiles)
  openr,1,filewithreducedimages
  readf,1,ReducedFiles
  close,1


  print,''
  print,'===ASTROMETRY==='
  print,''

;==Masking out extraordinary stars
  fakesky,ReducedFiles[0],calcpolout


;==taking the fakesky.fits image to get wcs and photometry
  imagewcs='fakesky.fits' 


;==Getting a list of the stars to calculate the WCS
  astrometrylist,imagewcs,calcpolout
  ;spawn,'rm '+calcpolout

;==Reading the (aproximate) ra and dec from image header
  imgwcs=readfits(imagewcs,hdrwcs)
  Ra_apx=sxpar(hdrwcs,'RA')
  Dec_apx=sxpar(hdrwcs,'DEC')

  coordstring=Ra_apx+' '+Dec_apx
  print,'Ra_apx , Dec_apx :',coordstring


  print,'Do you have Astrometry installed? (y/n):',astrome
  if (astrome eq 'y') then begin


;==Transform RA and Dec to degrees 
     stringad,coordstring,ra,dec
     imgsize=size(imgwcs)
     FOV=0.01* imgsize[2]       ;(arcmin) Approx size of an image 
     print,'Ra_apx , Dec_apx, FOV_apx :',ra,dec,FOV
     FOV_str = STRING(FOV)
     ra_str = STRING(ra)
     dec_str = STRING(dec)

     ;print,ra_str,dec_str,FOV_str


;==Using Astronmetry with fakesky.fit==
     spawn,'solve-field --ra '+ra_str+', --dec '+dec_str+', --radius '+FOV_str+' --overwrite fakesky.fits'
     spawn,'mv fakesky.new fakesky.new.fits'

  endif else begin
;==For no Astrometry.net===
     print,'Enter either the textXYsorted.list file or the fakesky.fits image in Astrometry.net'
     print,'to solve the astrometry'
     print,''
     wcsfile=''
     read,'Input the wcs.fits or the new-image.fits (output of Astrometry.net):',wcsfile

;==Put the wcs into the header
     hdr_wcs=headfits(wcsfile)
     extast,hdr_wcs,astr                   ;extract the astrometric information 
     putast,hdrwcs,astr                    ;put the astrometry into the hdr
     writefits,'fakesky.new.fits',imgwcs,hdrwcs ;write the image with the complete wcs info.
     print,'===Output fakesky.new.fits==='
  endelse

  wcsimage='fakesky.new.fits'
  image=ReducedFiles[0]

;==getting the name (removing the .fits)
  lengthname=strlen(ReducedFiles[0])
  namesdim=lengthname-5
  polfile=strmid(ReducedFiles[0],0,namesdim)+'.pol'


;==Catalogue with Ra dec and Mag===(solvepol second part)
;==Reading central position of wcsimage to get Ra & Dec and radio
  imagebg=readfits(wcsimage,hdr)
  sizex = n_elements(imagebg[*,0])
  sizey = n_elements(imagebg[0,*])

  xyad,hdr,sizex/2,sizey/2,Rac,Decc,/CELESTIAL
  xyad,hdr,0,0,Ra00,Dec00,/CELESTIAL
  xyad,hdr,sizex,0,Ra10,Dec10,/CELESTIAL
  xyad,hdr,0,sizey,Ra01,Dec01,/CELESTIAL


;==Radius of image==
  cosA = sin(Decc*!pi/180.)*sin(Dec00*!pi/180.) + cos(Decc*!pi/180.)*cos(Dec00*!pi/180.)*cos((Rac-Ra00)*!pi/180.)
  A = 180./!pi*acos(cosA)       ;degrees
  r_img = A * 60.0              ;radius of image in arcmin.


;==Geting catalogue and giving it structure
  info=QueryVizier('GSC2.3',[Rac,Decc],r_img,/ALLCOLUMNS)
  print,'Query Vizier site: France. Catalogue: GSCv2.3.'

  PRINT_STRUCT,info,labelarray,labelstr,catdata,FORM_FLOAT = ['f','11','6'],/STRING

;==Printing catalogue file
  catalog=subdir+'/'+'GSC2_3.txt'
  nlinesCat=N_ELEMENTS(catdata[0,*])
  close,1
  openw,1,catalog 

  ;label: VMAG , E_VMAG , RAJ2000 , DEJ2000
  printf,1,labelarray[15],' ',labelarray[16],' ',labelarray[4],' ',labelarray[5] 

  for i=0, nlinesCat-1 do begin
     printf,1,catdata[15,i],catdata[16,i],catdata[4,i],catdata[5,i]
  endfor
  close,1


  
;==Calibrating the magnitude==
  lengthname=strlen(ReducedFiles[0])
  namesdim=lengthname-5
  pairfile=strmid(ReducedFiles[0],0,namesdim)+'.pair'

  magnit,image,wcsimage,file=pairfile,cat=catalog

;==Plotting mag_measured vs mag_cat
  wcstest,subdir+'/'+'magnitall.out',catalog


;==Producing final table
  fintab,subdir+'/'+'magnitsigmaall.out',subdir+'/'+polfile


;==Plotvec === 
  finaltable=subdir+'/'+name+'fintab.out'
  readcol,finaltable,id,X,Y,Ra,Dec,ra_hrs,ra_min,ra_sec,dec_deg,dec_min,dec_sec,mag,dmag,p,dp,theta,dtheta,F_o_sigma,F_e_sigma,$
          format='I,D,D,D,D,i,i,f,i,i,f,D,D,D,D,D,D,D,D',COUNT=ngood,NLINES=nlines,/silent

  appa=90.                      ; in degrees
  print,'Over plot the ID of the stars on the polarisation map? (y/n):',overplotID

  IF (overplotID eq 'y') THEN plotvec,wcsimage,X,Y,p,theta,appa=appa,idstar=id $ 
  else plotvec,wcsimage,X,Y,p,theta,appa=appa


;==Ploting magnitude distribution & pol vs mag. 
  !P.MULTI = [0, 3, 1] 
  Device, RETAIN=2
  window,10,xsize=1050,ysize=300

  cgHistoplot,mag,xtitle='Magnitude',axiscolorname='White',color='White',charsize=3,binsize=0.5,BACKCOLORNAME='Black'
  
  plot,mag,p*100,PSYM=1,xtitle='Magnitude',ytitle='Polarisation (%)',charsize=3
  oploterror,mag,p*100,dmag,dp*100

  plot,mag,dp*100,PSYM=1,xtitle='Magnitude',ytitle='sigma Polarisation (%)',charsize=3
  !P.MULTI = 0 


;==Ploting magnitude distribution & pol vs mag to a eps file=====. 
  set_plot, 'ps'
  Device,filename=subdir+'/'+name+'magnitude.eps',XSIZE=20,YSIZE=20,/ENCAPSULATED
  cgHistoplot,mag,xtitle='Magnitude',axiscolorname='black',color='black',charsize=2,binsize=0.5,thick=5,CHARTHICK=5

  Device,filename=subdir+'/'+name+'magvspol.eps',XSIZE=20,YSIZE=20,/ENCAPSULATED
  plot,mag,p*100,PSYM=1,xtitle='Magnitude',ytitle='Polarisation (%)',charsize=2,xthick=5,ythick=5,CHARTHICK=5,thick=5
  oploterror,mag,p*100,dmag,dp*100,PSYM=1,thick=5

  Device,filename=subdir+'/'+name+'spolvsmag.eps',XSIZE=20,YSIZE=20,/ENCAPSULATED
  plot,mag,dp*100,PSYM=1,xtitle='Magnitude',ytitle='sigma Polarisation (%)',charsize=2,xthick=5,ythick=5,CHARTHICK=5,thick=5
  device,/close,ENCAPSULATED=0
  set_plot,'x'

;==find standard Star on the FOV and in the final catalogue==
  stdstars=solvepolDir + 'polstdstars.dat'
  findstdstar,subdir+'/'+pairfile,wcsimage,stdstars,RA,DEC

  print, '==== SOLVEPOL succesfully finished.  ==='

;;Ask for coffee
;coffee

  journal
end
