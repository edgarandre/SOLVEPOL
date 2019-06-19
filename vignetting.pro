pro vignetting,object

subdir='solvepol'

name=strmid(object,0,strlen(object)-5)

nlines = FILE_LINES(subdir+'/'+name+'.find')

image=readfits(object)
sizex = n_elements(image[*,0])
sizey = n_elements(image[0,*])


center=[sizex/2.,sizey/2.]
radius=sizex/2.-5.

header=strarr(13)

openr,1,subdir+'/'+name+'.find'
readf,1,header
close,1
readcol,subdir+'/'+name+'.find',F='I,D,D,D,D,D',starnum,xval,yval,flux,sharpness,roundness,skipline=13,count=ngood,/SILENT



print,'No. of detected sources: ',ngood

;Max number of permited sources (Int): 32767
if ngood gt 32767 then begin
   print,'The number of stars is greater than the maximum allowed (32767). Try a higher value for flux/sigma or set the FWHM to 10.'
   retall
endif

;;square which only includes points within border px from the edge
cuttype='square'
border=20
cutvals=where(yval(*) gt border and yval(*) lt sizey-border and xval(*) gt border and xval(*) lt sizex-border, numcutvals)

nvingetting_reject=ngood-numcutvals
print,'No. of sources rejected by VIGNETTING criteria: ',nvingetting_reject
print,'No. of detected sources after VIGNETTING: ',numcutvals



openw,1,subdir+'/'+name+'.cut'
printf,1,transpose(header[0:n_elements(header)-2])
printf,1,''
printf,1,' Program: VIGNETTING_CUT_VER2 '+systime()
printf,1,''
printf,1,' Cut performed is: ',cuttype
printf,1,' Cut paramters:', border,'< x < ',sizex-border ,'  ',border ,'< y < ', sizey-border
printf,1,''
printf,1,' No. of sources rejected by VIGNETTING criteria',nvingetting_reject
printf,1,''
printf,1,header[n_elements(header)-1]
format='(12x,i5,2f8.2,f9.1,2f9.2)'
for j=0,numcutvals-1 do begin
    i=cutvals(j)
    printf,1,f=format,starnum[i],xval[i],yval[i],flux[i],sharpness[i],roundness[i]
endfor
close,1

end
