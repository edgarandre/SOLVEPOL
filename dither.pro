
;+
; NAME:
;       DITHER
; PURPOSE:
;        combine images  
; EXPLANATION:
;      IMGSFILE is a text file with the names of the images to be
;      combined. NDITHER is the number of images per waveplate
;      position. The procedure find the shift between dithered images
;      and  combine them. The output images have  the
;      size of the first image per position of the waveplate. 
; CALLING SEQUENCE:
;       DITHER,imgsfile,ndither=ndither  
;
; INPUTS:
;       imgsfile - Name of text file containing the names of the
;                  images to combine in a column.                    
;       ndither - numer of images per waveplate position. 
; OPTIONAL INPUT KEYWORDS:
;       
;
; OUTPUTS:
;        One combined image per waveplate position. 
;
; OPTIONAL OUTPUT KEYWORDS:
;       
;
; EXAMPLES:
;       File 'images.list' contains a list of 16 images. The first 4
;       images were taken on the first position of the waveplate,
;       the consecutive 4 images were taken on the second position of
;       the waveplate and so on.
;
;       IDL>  DITHER,'images.list',dither=4 
;       Output --> imavgcombine1.fits
;                  imavgcombine2.fits 
;                  imavgcombine3.fits
;                  imavgcombine4.fits 
;
; RESTRICTIONS:
;       
;
; PROCEDURES CALLED
;       IMCOMB (included here)
;
; REVISION HISTORY:
;       
;-


pro imcomb,imagesnames,xshift=xshift,yshift=yshift,nwp=n

  noImages=n_elements(imagesnames) ;number of images
  print,'Images:',imagesnames
  print,'Number of images:',noImages
  print,'Waveplate position:',n

  refimg=readfits(imagesnames[0],hdr0,/SILENT) ;taking size of first image  
  sizex = n_elements(refimg[*,0])
  sizey = n_elements(refimg[0,*])
  int=fltarr(sizex,sizey,noImages)


  ;===reading the images (arrays)===
  for i=0,noImages-1 do begin
     int[*,*,i]=readfits(imagesnames[i],hdr,/silent)
  endfor


  ;;===taking the biggest size image===
  maxXshift=max(ABS(xshift))
  maxYshift=max(ABS(yshift))
  newint=fltarr(sizex+2*maxXshift,sizey+2*maxYshift,noImages,/NOZERO)
  imgcombine=fltarr(sizex+2*maxXshift,sizey+2*maxYshift)


 ;;===taking the size of the first image===
 ; newint=fltarr(sizex,sizey,noImages,/NOZERO)
 ; imgcombine=fltarr(sizex,sizey)


  ;===shifting images reference to the first image=== 
  for i=0,noImages-1 do begin
     for j=0,sizex-1 do begin
        for k=0,sizey-1 do begin         

;           newint[j+xshift[i],k+yshift[i],i]=int[j,k,i]
           newint[j+maxXshift+xshift[i],k+maxYshift+yshift[i],i]=int[j,k,i]

        endfor
     endfor
  endfor

  ;===adding images=== 
  for j=0,sizex+2*maxXshift-1 do begin
     for k=0,sizey+2*maxYshift-1 do begin           
;        imgcombine[j,k]= median(newint[j,k,*],/even)
        imgcombine[j,k]= avg(newint[j,k,*])
     endfor
  endfor  


;;==Name of file [0] to add as prefix===
  imagename=strmid(imagesnames[0],0,strlen(imagesnames[0])-5)

  NSTR=STRTRIM(n, 2)
;  hdr0=headfits(imagesnames[0],/silent) ;modify header indicating
;  reduction (note. headfits changes the bits/value of the pixels)
  sxaddpar,hdr0,'HISTORY','Images avg combined:'
  for i=0, noImages-1 do begin
     sxaddpar,hdr0,'HISTORY',imagesnames[i]
  endfor
  ;==writing and trimming to the size of first image
  writefits,imagename+'avgcomb_'+NSTR+'.fits',imgcombine[maxXshift:sizex-1,maxYshift:sizey-1],hdr0
  print,'Ouput image: '+imagename+'avgcomb_'+NSTR+'.fits'
  print,''


end


pro dither,imgsfile,dither=ndither

;  imgsfile='imgs.list'          ;text file with all images 
;  ndither=4                     ;number of images per waveplate position

;  if n_elements(ndither) eq 0 then ndither=1


  ;===reading text file with images===
  noimgs = FILE_LINES(imgsfile)      ;Number of files 
  imgsnames = strarr(fix(noimgs[0])) ;define number of files array
 
  
  close,1
  openr,1,imgsfile              ;Open and read text file
  readf,1,imgsnames             ;names
  close,1
  nwp=1
  for n=0,noimgs-1,ndither do begin

   
  imageA=readfits(imgsnames[n],hrdA,/SILENT) ;Image of reference (First image)
  sizeimageA=size(imageA)
  x_cntrd_shifts=fltarr(ndither)
  y_cntrd_shifts=fltarr(ndither)

  for i=n+1,n+ndither-1 do begin

     print,'Shift between '+imgsnames[n]+' and '+imgsnames[i]+' is:'
     imageB=readfits(imgsnames[i],hdrB,/SILENT) ;Shifted Images 
     sizeimageB=size(imageB)

     ;1st iteration
     correl_optimize,imageA[20:sizeimageA[1]-20,20:sizeimageA[2]-20],imageB[20:sizeimageB[1]-20,20:sizeimageB[2]-20],x_cntrd_shift,y_cntrd_shift,/NUMPIX,magnification=2

     ;2nd iteration
     ;correl_optimize,imageA[20:sizeimageA[1]-20,20:sizeimageA[2]-20],imageB[20:sizeimageB[1]-20,20:sizeimageB[2]-20],x_cntrd_shift,y_cntrd_shift,XOFF_INIT=x_cntrd_shift,YOFF_INIT=y_cntrd_shift,/NUMPIX,magnification=2

     
     ;x_cntrd_shift = 0
     ;y_cntrd_shift = 0
     
     print,'       x shift = ',x_cntrd_shift
     print,'       y shift = ',y_cntrd_shift
     
     x_cntrd_shifts[i-n]= x_cntrd_shift
     y_cntrd_shifts[i-n]= y_cntrd_shift

  endfor


  imcomb,imgsnames[n:n+ndither-1],xshift=x_cntrd_shifts,yshift=y_cntrd_shifts,nwp=nwp++

endfor


  print,'Creating text file with avg combined images: objectscombined.list'
  spawn,'ls *avgcomb* >  objectscombined.list'


end
