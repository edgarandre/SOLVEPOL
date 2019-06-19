pro wcstest,magnitall,catfile
;+
; NAME:
;
;      WCSTEST
;
; PURPOSE:
;
;      TESTS THE WCS BY MATCHING THE MEASURED POSITIONS AND THE
;      POSITIONS ON THE CATALOGUE   
;      
; EXPLANATION:
;        
;	Plots the measured magnitud and the magnitude from the given
;	catalog, and compare them with a linear line. 
;	
;
; CALLING SEQUENCE:
;
;     wcstest,magnitall,catfile
;
; INPUTS:
;
;       magnitall - output file of magnit.pro with ID,X,Y,Ra,Dec,Mag
;                   before p/sigma
;
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
;       wcstest,'/solvepol/magnitall.out','/solvepol/GSC2_3.txt'
;
;       
; PROCEDURES USED:
;       
; NOTES:
;       
; REVISON HISTORY:
;       
;-  	


  COMMON imgnames,name
  subdir='solvepol'
  !P.MULTI = 0 

;;--Read magnitall.out file from magnit.pro
  close,1
  openr,1,magnitall
  nolines= FILE_LINES(magnitall)

 ;read header (comment lines)------
  nHeader=2
  header=strarr(nHeader)
  readf,1,header

  mydata=dblarr(6,nolines-nHeader)
  readf,1,mydata
  close,1
;==note: RA=mydata[3,*],DEC=mydata[4,*]

;;--Read catalgue with Mag,Magerr,RA,Dec--
  readcol,catfile,F='D,D,D,D',MagCat,MagCatErr,RaCat,DecCat,COUNT=ngoodCat,NLINES=nlinesCat,/silent,COMMENT='N'


;;open a temporal file to write the matched magnitudes given Ra & Dec
  output='tmp_to_plot.txt'
  close,3
  openw,3,output
  sigmaRA=0.001                 ;uncertainty range to find the matching positions 
  sigmaDEC=0.001

;WHERE(FINITE(MagCatErr), COUNT = count).
  for i = 0,nolines-nHeader-1 do begin
     for j=0,ngoodCat-1 do begin
        if (mydata[3,i] ge (RaCat[j]-sigmaRA)) and (mydata[3,i] le (RaCat[j]+sigmaRA)) $
           and (mydata[4,i] ge (DecCat[j]-sigmaDEC)) and (mydata[4,i] le (DecCat[j]+sigmaDEC)) $
           and (FINITE(MagCatErr[j],/NaN) ne 1)  then begin

           printf,3,fix(mydata[0,i]),mydata[5,i],MagCat[j],MagCatErr[j]

        endif
     endfor
  endfor

  close,3

;;read the temporal file to produce the plot of mag-measured vs mag-catalog
  readcol,output,Format='I,D,D,D',IDmymag,mymag,catmag,catmagerr,COUNT=ngood,NLINES=nlines,/silent


  if N_elements(catmagerr) ne 0 then begin
     medianMagCatErr=median(catmagerr,/even)
     print,'Median of the Magnitude error in the GSC2.3: ',medianMagCatErr
  endif else begin
     medianMagCatErr=!Values.F_NAN
     print,'There are not matching stars in the GSC2.3 to estimate'
     print,' the Median of the Magnitude error in the catalogue'
     print,'Median of the Magnitude error in the GSC2.3: ',medianMagCatErr
  endelse

;;if there are no matching stars then fix magsigma=0
  if (ngood le 3) then begin
     print,''
     print,"There aren't enough stars matching the cataloge to plot Measured magnitude vs Catalog magnitude"
     print,"and calculate the magnitude error."
     print,''
     print,'Fixing the final Magnitude error to median mag error reported in the GSC2.3 = ',medianMagCatErr
     std_dev=medianMagCatErr
     
     ;std_dev=0.0
     ;print,'Final Magnitude error (Standard deviation)=',std_dev


     close,6
     openw,6,subdir+'/'+'magnitsigmaall.out'
     printf,6,header[0]
     printf,6,header[1] + '     MAGSIGMA'


     for i = 0,nolines-nHeader-1 do begin
        printf,6,mydata[0,i],mydata[1,i],mydata[2,i],mydata[3,i],mydata[4,i],mydata[5,i],std_dev,format='(I5, 2(F8.2," "),F9.5," ",F10.6," ",2(F9.5," "))'
     endfor
     
     close,6

     spawn,'rm tmp_to_plot.txt'

     return
  endif 


;=======FIRST INTERACTION TO ESTIMATE MAGERR===============
  catmagsortID=bsort(catmag,catmagsort)

;;===fit mymag catmag
  coeff=poly_fit(catmagsort,mymag[catmagsortID],1,yfit=yfitmags,sigma=sigma_coeff)

;;--Estimating dispersion std_deviation--
  deltamag=(coeff[0]+coeff[1]*catmag)-mymag
  npoints = N_ELEMENTS(deltamag)

  std_dev = SQRT(total((deltamag)*(deltamag))/npoints)
  print,'Magnitude standard deviation 1 =',std_dev


;;--Plot the magnitudes firs interaction--
  x=findgen(20)

;;===fitting a Gaussian to compare the FWHM with the std_dev;======================
  yhist=HISTOGRAM(deltamag,binsize=0.05,locations=xhist)
  if N_elements(yhist) gt 3 then begin
     Result = GAUSSFIT(xhist,yhist,Gauss_coeff,sigma=sigma,NTERMS=3)
     FWHM_gauss = 2*SQRT(2*ALOG(2))*Gauss_coeff[2]
     print,'Magnitude FWHM_gauss 1 =', FWHM_gauss
  endif

 ;;==plottig MeasuredMag vs CatMAg Excluding  magnitudes out of 3*std_dev======================
  output2='tmp_to_plot_2.txt'
  close,4
  openw,4,output2

  for i=0,npoints-1 do begin
     if (sqrt(deltamag[i]*deltamag[i]) LT 3*std_dev) then begin
        printf,4,IDmymag[i],mymag[i],catmag[i],catmagerr[i]
     endif
  endfor
  close,4

;;==read the temporal file to produce the plot of mag-measured vs mag-catalog
  readcol,output2,Format='I,D,D,D',IDmymagclean,mymagclean,catmagclean,catmagcleanerr,COUNT=ngoodclean,NLINES=nlinesclean,/silent


  if ngoodclean le 3 then begin
     print,''
     print,"There aren't enough stars matching the cataloge for the second interaction"
     print,"to plot Measured magnitude vs Catalog magnitude."
     print,'Final Magnitude error (Standard deviation 1)=',std_dev


     if std_dev le medianMagCatErr then begin
        print,'The Final Magnitude error is lower thant the median error reported in the GSC2.3: ',medianMagCatErr
        std_dev=medianMagCatErr
        print,'Fixing the Final Magnitude error to median mag error reported in the GSC2.3 = ',std_dev
     endif


     close,6
     openw,6,subdir+'/'+'magnitsigmaall.out'
     printf,6,header[0]
     printf,6,header[1] + '     MAGSIGMA'
     

     ;==keeping all mags fixing signama mag to the std_dev

     for i = 0,nolines-nHeader-1 do begin
        printf,6,mydata[0,i],mydata[1,i],mydata[2,i],mydata[3,i],mydata[4,i],mydata[5,i],std_dev,format='(I5, 2(F8.2," "),F9.5," ",F10.6," ",2(F9.5," "))'
     endfor

     close,6

     spawn,'rm tmp_to_plot_2.txt'

     return

  endif 
;=======SECOND INTERACTION TO ESTIMATE MAGERR===============
 
;;fit mymag catmag
  coeffclean=poly_fit(catmagclean,mymagclean,1,yfit=yfitmagsclean,sigma=sigma_coeffclean)

;;==Estimating dispersion
  deltamagclean=(coeffclean[0]+coeffclean[1]*catmagclean)-mymagclean

  npointsclean=N_ELEMENTS(deltamagclean)
  std_devclean = SQRT(total((deltamagclean)*(deltamagclean))/npointsclean)
  print,'Magnitude standard deviation 2 =',std_devclean


;;fitting a Gaussian to compare the FWHM with the std_dev;======================
  yhistclean=HISTOGRAM(deltamagclean,binsize=0.05,locations=xhistclean)
  if N_elements(yhistclean) gt 3 then begin
     Resultclean  = GAUSSFIT(xhistclean,yhistclean,Gauss_coeffclean,sigma=sigmaclean,NTERMS=3)
     FWHM_gaussclean = 2*SQRT(2*ALOG(2))*Gauss_coeffclean[2]
     print,'Magnitude FWHM_gauss 2 =', FWHM_gaussclean
  endif


;;Excluding  magnitudes out of 3*std_devclean======================
  output3='tmp_to_plot_3.txt'
  close,5
  openw,5,output3


  for i=0,npointsclean-1 do begin
     if (sqrt(deltamagclean[i]*deltamagclean[i]) LT 3*std_devclean) then begin
        magsclean2=[IDmymagclean[i],mymagclean[i],catmagclean[i],catmagcleanerr[i]]
        printf,5,magsclean2
     endif
  endfor

  close,5

;;read the temporal file to produce the plot of mag-measured vs mag-catalog
  readcol,output3,Format='I,D,D,D',IDmymagclean2,mymagclean2,catmagclean2,catmagcleanerr2,COUNT=ngoodclean2,NLINES=nlinesclean2,/silent



  if ngoodclean2 le 3 then begin
     print,''
     print,"There aren't enough stars matching the cataloge for the third interaction"
     print,"to plot Measured magnitude vs Catalog magnitude."
     std_dev=std_devclean
     print,'Final Magnitude error (Standard deviation 2)=',std_dev


    if std_dev le medianMagCatErr then begin
        print,'The Final Magnitude error is lower thant the median error reported in the GSC2.3: ',medianMagCatErr
        std_dev=medianMagCatErr
        print,'Fixing the Final Magnitude error to median mag error reported in the GSC2.3 = ',std_dev
     endif

     close,6
     openw,6,subdir+'/'+'magnitsigmaall.out'
     printf,6,header[0]
     printf,6,header[1] + '     MAGSIGMA'
     
     for i = 0,nolines-nHeader-1 do begin
        printf,6,mydata[0,i],mydata[1,i],mydata[2,i],mydata[3,i],mydata[4,i],mydata[5,i],std_dev,format='(I5, 2(F8.2," "),F9.5," ",F10.6," ",2(F9.5," "))'
     endfor

     close,6

     spawn,'rm tmp_to_plot_3.txt'

     return

  endif 


;;=======THIRD INTERACTION TO ESTIMATE MAGERR===============
  Device,Retain=2
  window,3,xsize=400,ysize=400
  plot,catmagclean2,mymagclean2,ytitle='Measured magnitude',xtitle='GSC2.3 magnitude',PSYM=1,CharSize=1.5,ystyle=16,background=0,COLOR=16777215


 ;;--fitting mymag catmag --
  coeffclean2=poly_fit(catmagclean2,mymagclean2,1,yfit=yfitmagsclean2,sigma=sigma_coeffclean2)

  oplot,x,coeffclean2[0]+coeffclean2[1]*x,LINESTYLE=0
  oplot,[0,5,10,15,20],[0,5,10,15,20],LINESTYLE=1 ;straight line to compare
  al_legend,['Equality','Linear fit'],LINESTYLE=[1,0],colors=['white','white'],outline_color ='white'

 ;;Estimating dispersion--
  deltamagclean2=(coeffclean2[0]+coeffclean2[1]*catmagclean2)-mymagclean2
  npointsclean2=N_ELEMENTS(deltamagclean2)

  std_devclean2 = SQRT(total((deltamagclean2)*(deltamagclean2))/npointsclean2)
  print,'Final Magnitude error (Standard deviation 3)=',std_devclean2

;;fitting a Gaussian to compare the FWHM with the std_dev;======================
  yhistclean2=HISTOGRAM(deltamagclean2,binsize=0.05,locations=xhistclean2)
  if N_elements(yhistclean2) gt 3 then begin
     Resultclean2  = GAUSSFIT(xhistclean2,yhistclean2,Gauss_coeffclean2,sigma=sigmaclean2,NTERMS=3)
     FWHM_gaussclean2 = 2*SQRT(2*ALOG(2))*Gauss_coeffclean2[2]
     print,'Magnitude FWHM_gauss 3 =', FWHM_gaussclean2
  endif


  if std_devclean2 le medianMagCatErr then begin
     print,'The Final Magnitude error is lower thant the median error reported in the GSC2.3: ',medianMagCatErr
     std_devclean2=medianMagCatErr
     print,'Fixing the Final Magnitude error to median mag error reported in the GSC2.3 = ',std_devclean2
  endif


;;===Printing the final interation==
  close,6
  openw,6,subdir+'/'+'magnitsigmaall.out'
  printf,6,header[0]
  printf,6,header[1] + '     MAGSIGMA'


;;keeping all mags within and out 3*std_devclean2 fixing signama mag to std_devclean2
;;--oploting mags within and out 3*std_devclean2 fixing magerr to std_devclean2-- 
  std_devclean2array = fltarr(N_ELEMENTS(mymag))
  std_devclean2array[*] = std_devclean2
  oploterror,catmag,mymag,catmagerr,std_devclean2array,psym=1,ERRCOLOR='White'
  for i = 0,nolines-nHeader-1 do begin

     printf,6,mydata[0,i],mydata[1,i],mydata[2,i],mydata[3,i],mydata[4,i],mydata[5,i],std_devclean2,format='(I5, 2(F8.2," "),F9.5," ",F10.6," ",2(F9.5," "))'

  endfor
  close,6


  spawn,'rm tmp_to_plot.txt'
  spawn,'rm tmp_to_plot_2.txt'
  spawn,'rm tmp_to_plot_3.txt'


;;;===== plotting to a eps file =====
  set_plot, 'ps'
  Device,filename=subdir+'/'+name+'magvscatmag.eps',XSIZE=20,YSIZE=20,/ENCAPSULATED
  plot,catmagclean2,mymagclean2,ytitle='Measured magnitude',xtitle='GSC2.3 magnitude',PSYM=1,YSTYLE=16,charsize=2,xthick=5,ythick=5,CHARTHICK=5,thick=5

  oplot,x,coeffclean2[0]+coeffclean2[1]*x,LINESTYLE=0,thick=5
  oplot,[0,5,10,15,20],[0,5,10,15,20],LINESTYLE=1,thick=5 ;straight line to compare

  al_legend,['Equality','Linear fit'],LINESTYLE=[1,0] ;,box=0,pos=[7,15]

  oploterror,catmag,mymag,catmagerr,std_devclean2array,psym=1,thick=5
  device,/close,ENCAPSULATED=0
  set_plot, 'x'


  print,'==== wcstest succesfully finished.  ==='

end
