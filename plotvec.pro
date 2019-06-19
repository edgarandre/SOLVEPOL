;+
; NAME:
;       PLOTVEC
;
; PURPOSE:
;      
;       OVERPLOT THE POLARISATION VECTORS OVER AN INTENSITY IMAGE
;
; EXPLANATION:
;
;       Oveplot p(theta) on a DSS image on the FOV of 'image' on the
;       positions X,Y. Image has to have wcs to get RA,DEC given
;       X,YRead FINTABs ASCII files to filter.
;
; CALLING SEQUENCE:
;        
;       plotvec,image,X,Y,p,theta,appa=appa,idstar=idstar
;
; INPUTS:
;
;       image - image of reference with the wcs information
;
;       X - vector with the X position in pixels 
;
;       Y - vector with the Y position in vectors
;
;       p - vector with the polarisation
;
;       theta - vector with the angle of the polarsiation
;
;       appa - apperture correction of the PA
;
;       idstar - vector with thte ID of the stars
;
; OPTIONAL INPUT KEYWORDS:
;       
;
; OUTPUTS:
;        
;       eps file with the polarisation vectors over the (r-band) SDSS image 
;
; OPTIONAL OUTPUT KEYWORDS:
;       
;
; EXAMPLES:
;       
;
;       
; RESTRICTIONS:
;       
;
; PROCEDURES CALLED
;       SUBCELL.PRO
;       PLOTIMAGE.PRO
;
; REVISION HISTORY:
;       
;-

pro plotvec,image,X,Y,p,theta,appa=appa,idstar=idstar

  COMMON imgnames,name
  subdir='solvepol'
;read image===
  image_ref=readfits(image,hdr)
  imSize=SIZE(image_ref)

  xyad,hdr,imSize[1]/2,imSize[2]/2,Rac,Decc,/CELESTIAL ;RA and DEC of centre
;xyad,hdr,0,imSize[2]/2,Ra0,Dec0,/CELESTIAL ;Ra and DEC of left centre side
  xyad,hdr,0,0,Ra00,Dec00,/CELESTIAL       ;RA and DEC of [0,0]
  xyad,hdr,imSize[1],0,Ra10,Dec10,/CELESTIAL ;RA and DEC of [x,0]
  xyad,hdr,0,imSize[2],Ra01,Dec01,/CELESTIAL ;RA and DEC of [0,y]


;cosA = sin(Dec0*!pi/180.)*sin(Decc*!pi/180.) + cos(Dec0*!pi/180.)*cos(Decc*!pi/180.)*cos((Ra0-Rac)*!pi/180.)
;A = 180./!pi*acos(cosA) ;Radius of image in degrees
;FOV = 2.* A * 60. ;size length of image in arcmin.

;;delta X of image==
  cosDX = sin(Dec00*!pi/180.)*sin(Dec10*!pi/180.) + cos(Dec00*!pi/180.)*cos(Dec10*!pi/180.)*cos(abs(Ra00-Ra10)*!pi/180.)
  DX=180./!pi*acos(cosDX)       ;degrees
  DX = DX * 60.0                ;delta X of image in arcmin.

;;delta Y of image==
  cosDY = sin(Dec00*!pi/180.)*sin(Dec01*!pi/180.) + cos(Dec00*!pi/180.)*cos(Dec01*!pi/180.)*cos(abs(Ra00-Ra01)*!pi/180.)
  DY=180./!pi*acos(cosDY)       ;degrees
  DY = DY * 60.0                ;delta Y of image in arcmin.

  FOV=max([DX,DY])              ;size of DSSimage = maximum size of the CCD image

  QueryDSS,[Rac,Decc],imagebg,Header,IMSIZE= FOV


;==set_plot,'x'
  !P.MULTI = 0 
  device, DECOMPOSED=0
  loadct, 39


;==defining window
  imagebgSize = SIZE(imagebg) 
  Window,5,title='DSS image (R-band)',xsize=700,ysize=600


  range=minmax(imagebg,subs,/NAN) 

  PLOTIMAGE,imagebg,range=range,/PRESERVE_ASPECT,XSTYLE=4,YSTYLE=4

  xyad,Header,0,0,Ra_0,Dec_0
  xyad,Header,imagebgSize[1],0,Ra_x,Dec_x
  xyad,Header,0,imagebgSize[2],Ra_y,Dec_y
  xyad,Header,imagebgSize[1],imagebgSize[2],Ra_xy,Dec_xy

  AXIS,XAXIS=0,XRANGE = [Ra_0,Ra_x],XSTYLE=1, XTITLE='RA (J2000)',/device
  AXIS,XAXIS=1,XRANGE = [Ra_y,Ra_xy],XSTYLE=1
  AXIS,YAXIS=0,YRANGE = [Dec_0,Dec_y],YSTYLE=1, YTITLE='DEC (J2000)',/device
  AXIS,YAXIS=1,YRANGE = [Dec_x,Dec_xy],YSTYLE=1

;=="weathervane" showing N-E orientation
  TVLCT, 255, 0, 0, 100         ;red
  arrows,Header,250,100,color=100,thick=2,arrowlen=4

  p_max=max(p,/NAN)
  barpixels =100                ;polarisation_max = 100 px
  length=barpixels/p_max              
  scale=1.0
  nlines=n_elements(X)

 ;==convert X and y to ra and dec using image header
  xyad,hdr, X, Y, RAs, DECs
 ;==convert ra and dec to X and Y using DSS image header
  adxy, Header, RAs, DECs, Xdss, Ydss

  TVLCT, 255, 255, 255, 101     ;white
  oplot,Xdss,Ydss,PSYM=3,thick=2  ;;optionally plot a cross on stars


  if  (n_elements(idstar) eq 0) then begin 

     for i=0, nlines-1 do begin

     ;projection in x and y of the p-vectors        
        dxdss=length*p[i]*cos(theta[i]*!dpi/180.0+appa*!dpi/180.0)
        dydss=length*p[i]*sin(theta[i]*!dpi/180.0+appa*!dpi/180.0)
      
     ;start and end of the p-vector
        strtxdss = Xdss - dxdss/2
        strtydss = Ydss - dydss/2
        endxdss = Xdss + dxdss/2
        endydss = Ydss + dydss/2

        oplot,scale*[strtxdss[i],endxdss[i]],scale*[strtydss[i],endydss[i]],thick=2,color=101 

        ;;test THE 90 DEGREES
        ;dxdsst=200*cos(90*!dpi/180.0+appa*!dpi/180.0)
        ;dydsst=200*sin(90*!dpi/180.0+appa*!dpi/180.0)
        ;strtxdsst = 300 - dxdsst/2
        ;strtydsst = 300 - dydsst/2
        ;endxdsst = 300 + dxdsst/2
        ;endydsst = 300 + dydsst/2
        ;oplot,scale*[strtxdsst,endxdsst],scale*[strtydsst,endydsst],thick=2,color=250  

     endfor

  endif else begin

     for i=0, nlines-1 do begin

     ;projection in x and y of the p-vectors        
        dxdss=length*p[i]*cos(theta[i]*!dpi/180.0+appa*!dpi/180.0)
        dydss=length*p[i]*sin(theta[i]*!dpi/180.0+appa*!dpi/180.0)
      
     ;start and end of the p-vector
        strtxdss = Xdss - dxdss/2
        strtydss = Ydss - dydss/2
        endxdss = Xdss + dxdss/2
        endydss = Ydss + dydss/2

        oplot,scale*[strtxdss[i],endxdss[i]],scale*[strtydss[i],endydss[i]],thick=2,color=101 

        XYOUTS,Xdss[i],Ydss[i],STRING(idstar[i],FORMAT='(i5)')

     endfor


  endelse 


;; Ploting the polarisation bar on top of the image
  oplot,[imagebgSize[1]*0.5,imagebgSize[1]*0.5+barpixels],[imagebgSize[2]*0.9,imagebgSize[2]*0.9],thick=3,color=255
  p_maxString = STRING(p_max*100, FORMAT='(F5.1)')
  xyouts, imagebgSize[1]*0.5, imagebgSize[2]*0.91, p_maxString+'%',CHARSIZE=2,CHARTHICK=2



;;=====plotting to a eps file =====
  set_plot, 'ps'
  loadct,0                      ;comment to colour

  Device,filename=subdir+'/'+name+'plotvec.eps', Decomposed=0,/portrait,/ENCAPSULATED;,XSIZE=20,YSIZE=20,/ENCAPSULATED
  plotimage,imagebg,range=reverse(range),/PRESERVE_ASPECT,XSTYLE=4,YSTYLE=4
;  plotimage,imagebg,range=range,/PRESERVE_ASPECT,XSTYLE=4,YSTYLE=4 ;colour


  TVLCT, 255, 255, 255, 101     ;white
  TVLCT, 0, 0, 0, 102           ;black 

  xyad,Header,0,0,Ra_0,Dec_0
  xyad,Header,imagebgSize[1],0,Ra_x,Dec_x
  xyad,Header,0,imagebgSize[2],Ra_y,Dec_y
  xyad,Header,imagebgSize[1],imagebgSize[2],Ra_xy,Dec_xy

  AXIS,XAXIS=0,XRANGE = [Ra_0,Ra_x],XSTYLE=1, XTITLE='RA (J2000)',xthick=5,CHARTHICK=5,/device
  AXIS,XAXIS=1,XRANGE = [Ra_y,Ra_xy],XSTYLE=1,xthick=5,CHARTHICK=5
  AXIS,YAXIS=0,YRANGE = [Dec_0,Dec_y],YSTYLE=1, YTITLE='DEC (J2000)',ythick=5,CHARTHICK=5,/device
  AXIS,YAXIS=1,YRANGE = [Dec_x,Dec_xy],YSTYLE=1 ,ythick=5,CHARTHICK=5

;;"weathervane" showing N-E orientation
  arrows,Header,0.4,0.2 ,color=102,thick=5,arrowlen=4,/NORMAL
;arrows,Header,0.4,0.2 ,thick=2,arrowlen=4,/NORMAL,color=100 ;colour

  oplot,Xdss,Ydss,PSYM=3,thick=5 ;;optionally plot a cross on stars

  for i=0, nlines-1 do begin

        ;projection in x and y of the p-vectors        
     dxdss=length*p[i]*cos(theta[i]*!dpi/180.0+appa*!dpi/180.0)
     dydss=length*p[i]*sin(theta[i]*!dpi/180.0+appa*!dpi/180.0)
      
        ;start and end of the p-vector
     strtxdss = Xdss - dxdss/2
     strtydss = Ydss - dydss/2
     endxdss = Xdss + dxdss/2
     endydss = Ydss + dydss/2

     oplot,scale*[strtxdss[i],endxdss[i]],scale*[strtydss[i],endydss[i]],thick=5,color=102
;        oplot,scale*[strtxdss[i],endxdss[i]],scale*[strtydss[i],endydss[i]],thick=2,color=101 ;colour

  endfor

;; Ploting the polarisation bar on top of the image
  oplot,[imagebgSize[1]*0.5,imagebgSize[1]*0.5+barpixels],[imagebgSize[2]*0.9,imagebgSize[2]*0.9],thick=5,color=102
  xyouts, imagebgSize[1]*0.5, imagebgSize[2]*0.91, p_maxString+'%',CHARSIZE=2,CHARTHICK=5, COLOR=102

;oplot,[imagebgSize[1]*0.5,imagebgSize[1]*0.5+100],[imagebgSize[2]*0.9,imagebgSize[2]*0.9],thick=3,color=101 ;colour
;xyouts, imagebgSize[1]*0.5, imagebgSize[2]*0.91, p_maxString+'%',CHARSIZE=2,CHARTHICK=3, COLOR=101 ;colour

  device,/close
  set_plot,'x'


  print, '==== plotvec succesfully finished.  ==='


end
