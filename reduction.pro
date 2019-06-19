pro reduction, objectfile, flats = flatfile, bias = biasfile
;+
; NAME:
;      REDUCTION
; PURPOSE:
;      DATA REDUCTION
; EXPLANATION:
;        setup path 
;	!PATH = Expand_Path('+/Users/Edgar/IDLWorkspace80/') + ':' + !PATH
;	
;       Reductionv2 applyied the standard procedure described in
;       Massey 1997 "A User's Gide to CCD Reduction with IRAF":
;   
;       'Bias median combined, overscan and sigma clipping filtered'
;       'Images bias and overscan corrected'
;       'Flats corrected by Bias, overscan, sigma clipping and normalized.'
;       images/flats_bias_combine and (1. / skymode)
;       'Triming reduced images:'
;       modify header indicating reduction
;
; CALLING SEQUENCE:
;     REDUCTION, objectfile, flats=, bias=
; INPUTS:
;     objectfile - text file listing image names
;
;     flats - text file listing flat names
;
;     bias - text file listing bias names
;
;     subdir - subdirectory where the reduced images will be
;              located. Sho
; OPTIONAL INPUTS:
;     
;
; OUTPUTS:
;    bias_zero.fits,  flat_bias_combine.fits, bias_flat_+objnames.fits
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


  biasim = readfits(biasfile)
  flatim = readfits(flatfile)

  ;==CCD detectors for direct camera on OPD==

  ;;CCD 105; 2048 x 2048

  ;;CCD 106; 1024 x 1024

  ;;CCD Andor iKon-L936-*; 2048 x 2048; overscan=[no overscan]
  xoverscan=0
  ovrscan1=0
  ovrscan2=0                 

  ;;CCD Andor iXon; 1024 x 1024

  ;;CCDTektronics; 1024 x 1024 ;overscan=[0:17,*]
  ;xoverscan=0
  ;ovrscan1=0
  ;ovrscan2=17                  

  ;;CCD 048; 800 x 1200; overscan= [1:770,1155:1200]
  ;;overscan along X-axis
  ;xoverscan=1
  ;ovrscan1=1155
  ;ovrscan2=1200
  
  ;;CCD 009

  ;;SITe 




  if N_ELEMENTS(flatim) eq 1 or N_ELEMENTS(biasim) eq 1 then begin ;Case the flat and bias are textfiles 

;===reading Bias text file===
     noBias = FILE_LINES(biasfile) ;Number of files 
     noBias=fix(noBias[0])
     biasnames = strarr(noBias) ;define number of files array

     openr,1,biasfile           ;Open and read text file
     readf,1,biasnames
     close,1
  
     refimgb = readfits(biasnames[0],hdrb) ;reference image 
     sizexb = n_elements(refimgb[*,0])
     sizeyb = n_elements(refimgb[0,*])
     zerocombine = fltarr(sizexb,sizeyb)
     intb = fltarr(sizexb,sizeyb,noBias)
     polyarraybias = fltarr(sizexb,sizeyb)
     axis_yb = FINDGEN(sizeyb)
     axis_xb = FINDGEN(sizexb)


     for i=0,noBias-1 do begin	;reading the arrays (bias)
        intb[*,*,i]=readfits(biasnames[i],/silent)
     endfor


;;== combining the bias==
     for j=0,sizexb-1 do begin
        for k=0,sizeyb-1 do begin
;           zerocombine[j,k] = median(intb[j,k,*])
           zerocombine[j,k] = avg(intb[j,k,*])
        endfor
     endfor

;;==average over all columns of the overscan region==
     if xoverscan eq 1 then begin

        bias_zero_ovrscan=fltarr(sizexb) 

        for i=0,sizexb-1 do begin
           for k=ovrscan1,ovrscan2-1 do begin
              bias_zero_ovrscan[i]=avg(zerocombine[i,k]) 
           endfor
        endfor
     

     endif else begin
     
        bias_zero_ovrscan=fltarr(sizeyb)


        for k=0,sizeyb-1 do begin
           for i=ovrscan1,ovrscan2-1 do begin
              bias_zero_ovrscan[k]=avg(zerocombine[i,k]) 
           endfor
        endfor

     endelse


;;==Fitting the overscan region of the bias==

     if xoverscan eq 1 then begin

        polyfitbias=poly_fit(axis_xb,bias_zero_ovrscan,2,yfit=yfitbias)
        for k=0,sizeyb-1 do begin
           for j=0,sizexb-1 do begin
              polyarraybias[j,k] = yfitbias[j]
           endfor
        endfor

     endif else begin

        polyfitbias=poly_fit(axis_yb,bias_zero_ovrscan,2,yfit=yfitbias)
        for j=0,sizexb-1 do begin
           for k=0,sizeyb-1 do begin
              polyarraybias[j,k] = yfitbias[k]
           endfor
        endfor

     endelse

;;==clean deviant pixels (sigma clipping)==

     print,''
     print,'Bias avg combined and sigma filtered' 
     zerocombine = sigma_filter(zerocombine,300,N_SIGMA=2.5,/MONITOR,/ALL_PIXELS)
     bias_zero = zerocombine - polyarraybias

     print,''
     print,'Creating bias_zero.fits - master bias -'
     print,''

     writefits,subdir+'/'+'bias_zero.fits',bias_zero, hdrb
     
     
;==Correcting flats by bias and overscan =============== 
     noFlats = FILE_LINES(flatfile) ;Number of flats 
     noFlats = fix(noFlats[0])
     flatnames=strarr(noFlats)  ;define array with number of flats
     
     openr,1,flatfile           ;Open text file with object names
     readf,1,flatnames
     close,1

     refimgf = readfits(flatnames[0],hdrf) ;reference Flat 
     sizexf = n_elements(refimgf[*,0])
     sizeyf = n_elements(refimgf[0,*])
     intf = fltarr(sizexf,sizeyf,noFlats)
     medianflat=fltarr(sizexf,sizeyf)
     flats_bias = fltarr(sizexf,sizeyf)


;==fit polynomio to overscanregion==
     axis_yf=FINDGEN(sizeyf)
     axis_xf=FINDGEN(sizexf)
;     yfitf=fltarr(noFlats)
     polyarrayf=fltarr(sizexf,sizeyf,noFlats)
     Flat_bias=fltarr(sizexf,sizeyf,noFlats)


;==Reading the flats===============
     for i=0,noFlats-1 do begin	

        intf[*,*,i]=readfits(flatnames[i],/silent) ;reading the Flats 

;==Correcting flat images by overscan =============== 

        if xoverscan eq 1 then begin

           flats_ovrscan=fltarr(sizexf,noFlats)

           for l=0,sizexf-1 do begin
              for k=ovrscan1,ovrscan2-1 do begin
                 flats_ovrscan[l,i]=avg(intf[l,k,i]) ;average over the overscan region (Y-scrunch)
              endfor
           endfor
 
           polyfitf=poly_fit(axis_xf,flats_ovrscan[*,i],2,yfit=yfitf) 

        endif else begin
           
           flats_ovrscan=fltarr(sizeyf,noFlats)

           for k=0,sizeyf-1 do begin
              for l=ovrscan1,ovrscan2-1 do begin
                 flats_ovrscan[k,i]=avg(intf[l,k,i]) ;average over the overscan region (X-scrunch)
              endfor
           endfor

           polyfitf=poly_fit(axis_yf,flats_ovrscan[*,i],2,yfit=yfitf) 

        endelse


;==Create array with resulted polynomia==
        if xoverscan eq 1 then begin

           for k=0,sizeyf-1 do begin
              for j=0,sizexf-1 do begin
                 polyarrayf[j,k,i]=yfitf[j]
              endfor
           endfor

        endif else begin
           
           for j=0,sizexf-1 do begin
              for k=0,sizeyf-1 do begin
                 polyarrayf[j,k,i]=yfitf[k]
              endfor
           endfor
           
        endelse

        Flat_bias[*,*,i] = intf[*,*,i]  - bias_zero  - polyarrayf[*,*,i]

     endfor


;==combining the flats===
     for j=0,sizexf-1 do begin
        for k=0,sizeyf-1 do begin
;           medianflat[j,k] = median(Flat_bias[j,k,*])
           medianflat[j,k] = avg(Flat_bias[j,k,*])
        endfor
     endfor

     medianflat=sigma_filter(medianflat, 300,N_SIGMA=2.5,/MONITOR,/ALL_PIXELS) ;apply sigma cliping
     print,''
     print,'Flats avg combined and sigma filtered'
     
     medianValFlat = median(medianflat) ;estimating median of flat field to normalize.
     print,'Median value of combined Flat field ',medianValFlat
     flats_bias = medianflat    ;/ medianValFlat  ;normalizing the mediancombined flat field


     print,'Flats corrected by Bias, overscan, median combined and sigma clipping.'
     print,''
     print,'Creating flat_bias_combine.fits - master flat -'
     print,''
     writefits,subdir + '/' + 'flat_bias_combine.fits', flats_bias, hdrf


  endif else begin              ;Case the flat and bias are images ============== ============== 

     bias_zero = biasim
     flats_bias = flatim

     medianValFlat = median(flats_bias) ;estimating median of flat field to normalize.
     print,'Median value of the Flat field ',medianValFlat

  endelse


;===Correcting objects by bias and overscan ==============
  noImages = FILE_LINES(objectfile) ;Number of images 
  noImages = fix(noImages[0])
  objnames=strarr(noImages)     ;define array with number of images

  openr,1,objectfile            ;Open text file with object names
  readf,1,objnames
  close,1

  refimg=readfits(objnames[0])	;reference image 
  sizex = n_elements(refimg[*,0])
  sizey = n_elements(refimg[0,*])
  int=fltarr(sizex,sizey,noImages)
  object_bias=fltarr(sizex,sizey,noImages)
  polyarray=fltarr(sizex,sizey,noImages)

;==fit polynomio to overscanregion==
  axis_y=FINDGEN(sizey)
  axis_x=FINDGEN(sizex)
;  yfit=fltarr(noImages)


  for i=0,noImages-1 do begin                    ;reading the arrays (images)
     
     int[*,*,i]=readfits(objnames[i],/silent)

     if xoverscan eq 1 then begin

        object_ovrscan=fltarr(sizex,noImages)

        for j=0,sizex-1 do begin
           for k=ovrscan1,ovrscan2-1 do begin
              object_ovrscan[j,i]=avg(int[j,k,i]) ;average over all columns of the overscan region (X-crunch)
           endfor
        endfor

        polyfit=poly_fit(axis_x,object_ovrscan[*,i],2,yfit=yfit)
        

     endif else begin
        
        object_ovrscan=fltarr(sizey,noImages)

        for k=0,sizey-1 do begin
           for j=ovrscan1,ovrscan2-1 do begin
              object_ovrscan[k,i]=avg(int[j,k,i]) ;average over all rows of the overscan region (Y-crunch)
           endfor
        endfor

        polyfit=poly_fit(axis_y,object_ovrscan[*,i],2,yfit=yfit)

     endelse

;;==Ploting the fit of the overscanregion== 
     window,1,xsize=700,ysize=300

     if xoverscan eq 1 then begin
        
        plot,axis_x,yfit,xtitle='overscan region',ytitle='Poly_fit',/YNOZERO,title='Fit of the overscanregion'

     endif else begin

        plot,axis_y,yfit,xtitle='overscan region',ytitle='Poly_fit',/YNOZERO,title='Fit of the overscanregion'

     endelse

;;create array with resulted polynomia
     if xoverscan eq 1 then begin

        for k=0,sizey-1 do begin
           for j=0,sizex-1 do begin
              polyarray[j,k,i]=yfit[j]
           endfor
        endfor

     endif else begin

        for j=0,sizex-1 do begin
           for k=0,sizey-1 do begin
              polyarray[j,k,i]=yfit[k]
           endfor
        endfor

     endelse

     object_bias[*,*,i] = int[*,*,i]  - bias_zero - polyarray[*,*,i] 

  endfor

  print,''
  print,'Images corrected by Bias and overscan'
  print,''


;==images/flats_bias_combine normalized==============

  object_bias_flat=fltarr(sizex,sizey,noImages)

  for i=0,noImages-1 do begin
     sky,object_bias[*,*,i],skymode,skysig

     object_bias_flat[*,*,i] = (object_bias[*,*,i]/skymode) / (flats_bias / medianValFlat)

     object_bias_flat[*,*,i] =  object_bias_flat[*,*,i] * skymode ; recuperating the sky value
  endfor


  print,''
  print,'Images corrected by Flats'
  print,''

  object_bias_flat_files='bias_flat_'+objnames

  print,'Triming and modifying header indicating reduction:'

  for i=0,noImages-1 do begin

     readhdr=readfits(objnames[i],hdr,/silent)
     sxaddpar,hdr,'Solvepol','Images processed by reduction.pro'

    ;writefits,object_bias_flat_files[i],object_bias_flat[,,i],hdr ;CCD 105; 2048 x 2048

    ;writefits,object_bias_flat_files[i],object_bias_flat[,,i],hdr ;CCD 106; 1024 x 1024

    writefits,object_bias_flat_files[i],object_bias_flat[*,*,i],hdr ;CCD Andor iKon-L936-*; 2048 x 2048

    ;writefits,object_bias_flat_files[i],object_bias_flat[,,i],hdr ;CCD Andor iXon; 1024 x 1024

    ; writefits,object_bias_flat_files[i],object_bias_flat[18:1022,0:1022,i],hdr ;CCDTektronics; 1024 x 1024
  
    ;writefits,object_bias_flat_files[i],object_bias_flat[*,0:1150,i],hdr ;CCD 048; 800 x 1200

    ;writefits,object_bias_flat_files[i],object_bias_flat[,,i],hdr ;CCD 009

    ;writefits,object_bias_flat_files[i],object_bias_flat[,,i],hdr ;SITe

    ;writefits,object_bias_flat_files[i],object_bias_flat[18:1040,0:1022,i],hdr ;Detector=WI101


     print,'Reduced image: ',object_bias_flat_files[i]
  endfor



;;==Creating file text with the reduced images===
  print,''
  print,'Text file with reduced names: objectsreduced.list'
  print,''
  openw,1,subdir + '/' + 'objectsreduced.list'

  for i=0,noImages-1 do begin
     printf,1,object_bias_flat_files[i]
  endfor
  close,1



;--tests of 2nd flat-field correction
;img=readfits('bias_flat_hd23512_L0_0medcomb_1.fits')

;ims=SMOOTH(img,5,/EDGE_TRUNCATE)
;ims=SMOOTH(img,100,/EDGE_TRUNCATE)

;ims = MEDIAN(img, 50)

;ims=sigma_filter(ims, 300,N_SIGMA=2.5,/MONITOR,/ALL_PIXELS)
;ims=sigma_filter(ims,200,N_SIGMA=2.,/MONITOR,/ALL_PIXELS)
;ims=sigma_filter(ims,100,N_SIGMA=1.,/MONITOR,/ALL_PIXELS)

;writefits,'ims.fits',ims

end
