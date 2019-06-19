pro fintab,magfile,polfile
;+
; NAME:
;
;      FINTAB
;
; PURPOSE:
;
;      CREATES THE FINAL CATALOGUE OR FINTAB.OUT 
;      
; EXPLANATION:
;        
;	creates a file with ' ID	      X	    Y	RA(J2000)
;	DEC(J2000)   RA(1950)       DEC(1950)       MAG     MAGSIGMA
;	P     PSIGMA	  THETA     TSIGMA        Flux_o/    Flux_e/'
;       using .par and magnit.out outputs files.
;	
;
; CALLING SEQUENCE:
;
;     fintab,magfile,polfile
;
; INPUTS:
;       
;       magfile - output file of magnit.pro with ID,X,Y,Ra,Dec,Mag
;
;
;       polfile - output of calcpol.pro with IDpol,Xpol,Ypol and
;                 q,dq,u,du,p,sigma,dp,theta,dtheta.
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
;      fintab,'solvepol/magnitsigmaall.out','solvepol/hd110cv0001.pol'
;
;       
; PROCEDURES USED:
;       
; NOTES:
;       
; REVISON HISTORY:
;
; Prints the J2000 RA & Dec only - the 1950 Ra & Dec have been
;                                  removed. 4 March 2019. EdgarARA.
;
;-  	

  COMMON sig,skysig                   ;from findsep.pro
  COMMON flux,starnu,fluxes,ngoodpair ;from .pair
  COMMON imgnames,name
  subdir='solvepol'

;===pol and theta === 
  readcol,polfile,format='I,D,D',IDpol,Xpol,Ypol,COUNT=ngoodpol,NLINES=nlinespol,/silent
  readcol,polfile,format='D,D,D,D,D,D,D,D,D',q,dq,u,du,p,sigma,dp,theta,dtheta,COUNT=ngood,NLINES=nlines,/silent


;magnit
  readcol,magfile,F='I,D,D,D,D,D,D',ID,X,Y,Ra,Dec,Mag,Magsigma,COUNT=ngoodMag,NLINES=nlinesMag,/silent 


;printing 
  close,1
  openw,1,subdir+'/'+name+'fintab.out'
   printf,1,' ID     X	   Y	RA(J2000.0) DEC(J2000.0) RA(J2000.0)    DEC(J2000.0)     MAG     MAGSIGMA   P     PSIGMA  THETA     TSIGMA     Flux_o/   Flux_e/'
  printf,1,'        (px)	  (px)	(degrees)   (degrees)   (HR,MIN,SEC)   (DEC,MIN,SEC)    (mag)     (mag)    (f)     (f)	  (degrees) (degrees)  sigma     sigma'
  for j=0,ngoodpol/3-1 do begin
     for i=0, ngoodMag-1 do begin 
        if (IDpol[j*3] eq ID[i]) then begin
           for k=0,ngoodpair/2-1 do begin
              if (IDpol[j*3] eq starnu[k*2]) then begin
                 printf,1,format='(I5," ", 2(F8.2," "),F9.5," ",F10.6," ",A," ",2(F8.3," "),2(F8.6," "),2(F8.4," "),2(f9.1," "))',$
                        ID[i],X[i],Y[i],RA[i],DEC[i],adstring([RA[i],DEC[i]],precision=4),Mag[i],magsigma[i],p[j],sigma[j],theta[j],dtheta[j],fluxes[k*2]/skysig,fluxes[k*2+1]/skysig
              endif
           endfor
        endif
     endfor
  endfor
  close,1


  print,'Fintab output: /'+subdir+'/'+name+'fintab.out'

  readcol,subdir+'/'+name+'fintab.out',format='I,D,D,D,D,i,i,f,i,i,f,D,D,D,D,D,D,D,D',idf,Xf,Yf,Raf,Decf,ra_hrsf,ra_minf,ra_secf,dec_degf,dec_minf,dec_secf,mag,dmagf,pf,dpf,thetaf,dthetaf,F_o_sigmaf,F_e_sigmaf,COUNT=ngoodf,NLINES=nlinesf,/silent



;;;plot on screen==========
  set_plot, 'x'
  !P.MULTI = [0, 2, 2] 
  Device, RETAIN=2
  window,4,XSIZE=700, YSIZE=600

;;plotting Q vs U =======
  plot,q,u,xtitle='Q',ytitle='U',PSYM=1,background=0,charsize=1.5
  oploterror,q,u,dq,du,PSYM=1,color='White',ERRSTYLE=1

;;Plotting pvectors======
  plot,Xf,Yf,xtitle='X [px]',ytitle='Y [px]',PSYM=3,xstyle=1,ystyle=1,background=0,charsize=1.5,COLOR=16777215

  p_max=max(pf,/NAN)
  length=100./p_max             ;polarisation_max = 100 px

  scale=1.0
  appa=0.0
  for i=0, ngoodf-1 do begin

        ;projection in x and y of the p-vectors        
     dx=length*pf[i]*cos(thetaf[i]*!dpi/180.0+appa*!dpi/180.0)
     dy=length*pf[i]*sin(thetaf[i]*!dpi/180.0+appa*!dpi/180.0)
      
        ;start and end of the p-vector
     strtx = Xf - dx/2
     strty = Yf - dy/2
     endx = Xf + dx/2
     endy = Yf + dy/2

     oplot,scale*[strtx[i],endx[i]],scale*[strty[i],endy[i]],thick=1,COLOR=16777215
       
  endfor

;;plotting histogram of Theta=======
  cgHistoplot,thetaf,xtitle='theta',axiscolorname='White',color='White',charsize=1.5,binsize=5.,BACKCOLORNAME='Black'


;;plotting histogram of polarisation======
  cgHistoplot,pf,xtitle='polarisation',axiscolorname='White',color='White',charsize=1.5,binsize=0.005,BACKCOLORNAME='Black'


  !P.MULTI = 0 
  print,'Number of Stars (Pairs) = ',ngoodf


;;;========== plotting  to eps ==========
  set_plot,'ps'

;;plotting Q vs U =======
  Device,filename= subdir+'/'+name+'QvsU.eps', XSIZE=20,YSIZE=20,/ENCAPSULATED
  plot,q,u,xtitle='Q',ytitle='U',PSYM=1 ,charsize=2,xthick=5,ythick=5,CHARTHICK=5,thick=5
  oploterror,q,u,dq,du,PSYM=1,ERRSTYLE=1,thick=5

;;Plotting pvectors======
  Device,filename=subdir+'/'+ name+'pvectors.eps', XSIZE=20,YSIZE=20,/ENCAPSULATED
  plot,Xf,Yf,xtitle='X [px]',ytitle='Y [px]',PSYM=3,xstyle=1,ystyle=1,charsize=2,xthick=5,ythick=5,CHARTHICK=5

  for i=0, ngoodf-1 do begin

        ;projection in x and y of the p-vectors        
     dx=length*pf[i]*cos(thetaf[i]*!dpi/180.0+appa*!dpi/180.0)
     dy=length*pf[i]*sin(thetaf[i]*!dpi/180.0+appa*!dpi/180.0)
     
        ;start and end of the p-vector
     strtx = Xf - dx/2
     strty = Yf - dy/2
     endx = Xf + dx/2
     endy = Yf + dy/2

     oplot,scale*[strtx[i],endx[i]],scale*[strty[i],endy[i]],thick=5
  endfor

;;plotting histogram of Theta=======
  Device,filename= subdir+'/'+name+'theta.eps', XSIZE=20,YSIZE=20,/ENCAPSULATED
  cgHistoplot,thetaf,xtitle='theta',axiscolorname='Black',color='Black',charsize=2,binsize=5.,thick=5,CHARTHICK=5

;;plotting histogram of polarisation======
  Device,filename= subdir+'/'+name+'polarisation.eps', XSIZE=20,YSIZE=20,/ENCAPSULATED
  cgHistoplot,pf,xtitle='polarisation',axiscolorname='Black',color='Black',charsize=2,binsize=0.005,thick=5,CHARTHICK=5

  device,/close,ENCAPSULATED=0
  set_plot, 'x'

  print, '==== fintab succesfully finished.  ==='

end
