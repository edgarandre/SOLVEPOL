pro runphot,filename,phpadu,apr,skyrad,badpix,fwhm
  subdir='solvepol'

  img=readfits(filename+'.fits',himg)

  readcol,subdir+'/'+filename+'.pair',F='I,D,D,D,D,D',starnum,xval,yval,flux,sharpness,roundness,COUNT=ngood,NLINES=nlines,/silent
	
  header=strarr(nlines-ngood-ngood/2)

  openr,1,subdir+'/'+filename+'.pair'
  readf,1,header
  close,1

  fwhm=fwhm
  gcntrd, img, xval, yval, xcen, ycen, fwhm,  /SILENT
  
  aper,img,xcen,ycen,flux,eflux,sky,skyerr,phpadu,apr,skyrad,badpix,/flux,/silent
;  aper,img,xval,yval,flux,eflux,sky,skyerr,phpadu,apr,skyrad,badpix,/flux,/silent


  
  openw,1,subdir + '/' +filename+'.phot'
  printf,1,transpose(header[0:n_elements(header)-2])
  printf,1,''
  printf,1,' Program: RUNPHOT '+systime()
  printf,1,''
  for i=0,n_elements(apr)-1 do printf,1,f='(a,i2,a,f4.1)','Radius of aperture ',i,' = ',apr[i]
  printf,1,f='(/a,f4.1)','Inner radius for sky annulus = ',skyrad[0]
  printf,1,f='(a,f4.1)', 'Outer radius for sky annulus = ',skyrad[1]
  printf,1,''
  printf,1,header[n_elements(header)-1]+' Sky    SkyErr    Flux     FluxErr'
	
  format='(12X,I5,2F8.2,F9.1,2F9.2,'+strtrim(n_elements(apr)*2+2,2)+'(F12.2))'
	
  for i=0,ngood/2-1 do begin
     printf,1,'PAIR',i+1
     j=i*2
     k=i*2+1

;     printf,1,f=format,starnum[j],xval[j],yval[j],flux[j],sharpness[j],roundness[j],sky[j],skyerr[j],flux[*,j],eflux[*,j]
;     printf,1,f=format,starnum[k],xval[k],yval[k],flux[k],sharpness[k],roundness[k],sky[k],skyerr[k],flux[*,k],eflux[*,k]

     printf,1,f=format,starnum[j],xcen[j],ycen[j],flux[j],sharpness[j],roundness[j],sky[j],skyerr[j],flux[*,j],eflux[*,j]
     printf,1,f=format,starnum[k],xcen[k],ycen[k],flux[k],sharpness[k],roundness[k],sky[k],skyerr[k],flux[*,k],eflux[*,k]
  endfor
  close,1
  
  return
end 
