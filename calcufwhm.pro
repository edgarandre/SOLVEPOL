function calculatethefwhm,wrefimg

  gauss_result=gauss2dfit(wrefimg,gauss_params)
  fwhm=sqrt(2.*alog(2.))*(gauss_params(2)+gauss_params(3))

  return,fwhm
  
end

pro calcufwhm,image,fwhm

  refimg=readfits(image)        ;read image
  srefimg=size(refimg)          ;size image
  windowsize=50                 ;side of square around stars to calculate fwhm

  trimrefimg=refimg[windowsize:srefimg[1]-windowsize,windowsize:srefimg[2]-windowsize] ;trimming borders

;  maxval= MAX(trimrefimg, location)


;;--Approach to find stars--
  sharplim=[0.2,1.0]            ;DEFAULT
  roundlim=[-1.0,1.0]           ;DEFAULT
  sky,refimg,skymode,skysig
  hmin =10. * skysig ;V16.1 is  3. * skysig
  fwhm=5.                       ;first guess of the fwhm

  find, trimrefimg, x, y, flux, sharp, roundness, hmin, fwhm, roundlim, sharplim


  if N_elements(flux) eq 0 then begin ;if there are no detected sources then fix to 5 the fwhm
     fwhm = 5.
     print,'The fixed FWHM is:', fwhm
     return
  endif

  fluxindx=where(flux gt 0)     ;avoid negative flux (valleys)
  flux=flux[fluxindx]
  sortmax=sort(-flux)
  nstars=N_elements(flux)       
  if nstars gt 30  then nstars = 30 ;number of stars to calculate fwhm up to 30

  sortmax=sortmax[0:nstars-1]
  x=x+windowsize
  y=y+windowsize
;;--

  nfwhm = 0

  nwindows=n_elements(sortmax)
  allfwhm=fltarr(nwindows)

  for i=0,nwindows-1 do begin   ;calculating the fwhm of the detected sources
 
     wrefimg=refimg[x[sortmax[i]]-windowsize:x[sortmax[i]]+windowsize-1,y[sortmax[i]]-windowsize:y[sortmax[i]]+windowsize-1]

     fwhm_i=calculatethefwhm(wrefimg)
     
     if fwhm_i gt 0 then begin
        print,format='(%"The fwhm at [%d, %d] is %f")',x[sortmax[i]],y[sortmax[i]],fwhm_i
     
        nfwhm++
        allfwhm[i]=fwhm_i
     endif
endfor


;;--the next method takes the 50% from the max peak--
; good = WHERE(trimrefimg ge maxval*0.5)
;  ind = ARRAY_INDICES(trimrefimg, good)
  
 ; ind=ind+50

;  nwindows=n_elements(ind[0,*])
;  allfwhm=fltarr(nwindows)

;  for i=0,nwindows-1 do begin

;     wrefimg=refimg[ind[0,i]-windowsize:ind[0,i]+windowsize-1,ind[1,i]-windowsize:ind[1,i]+windowsize-1]

;     maxval2=max(wrefimg)

;     if maxval2 eq refimg[ind[0,i],ind[1,i]] then begin

;        fwhm_i=calculatethefwhm(wrefimg)

;        print,format='(%"The fwhm at [%d, %d] is %f")',ind[0,i],ind[1,i],fwhm_i

;        nfwhm++
;        allfwhm[i]=fwhm_i

;     endif
;  endfor



  print,'Number of stars to calculate the fwhm',nfwhm
  allfwhmindex = where(allfwhm ne 0)
  allfwhm=allfwhm[allfwhmindex]

;  window,xsize=400,ysize=400
;  y_hfwhm = HISTOGRAM(allfwhm,binsize=1,locations=x_hfwhm)
;  plot,x_hfwhm,y_hfwhm,psym=10,xtitle='FWHM',ytitle='Number density'

   ;--FWHM mode--
;  maxfreq = Max(y_hfwhm)
;  fwhm = Where(y_hfwhm eq maxfreq, count) + min(allfwhm)
;  print,'The mode FWHM is:', fwhm

  ;--FWHM average--
  fwhm = avg(allfwhm)
  print,'The average FWHM is:', fwhm

  ;--FWHM median--
;  fwhm = median(allfwhm)
;  print,'The median FWHM is:', fwhm


;  if allfwhmindex eq -1 then fwhm = 10 ;if there are not stars to measure the fwhm (unlikely), the fix to 10 the fwhm.

end
