pro calcpol,filenames,nwp,rdn,deltatheta

  COMMON Polsigma,p_over_sigma,Chisqrt_lim,chisqrt_int
  !EXCEPT=0
  subdir ='solvepol'

  JUMP1: PRINT, 'Calculating with the new Chisqrt'

  ;read,'Give me pol/sigma_pol value : ', p_over_sigma
  print,'pol/sigma_pol value: ',p_over_sigma

;Test to make sure the number of files equals the number of waveplate positions
  if n_elements(filenames) ne nwp then begin
     print,'ERROR: The number of files (',strtrim(n_elements(filenames),2),') is not the same as the number of waveplate positions (',strtrim(nwp,2),')'
     print,'...Quitting.'
     retall
  endif


  napr=0                        ;number of apertures
  tmp=''                        ;used to read in lines from <filename>.phot
  openr,1,subdir+'/'+filenames[0]+'.phot'
  readf,1,tmp
  nhlines=1                     ;number of header lines
  print,tmp
  while ~strcmp(tmp,'             STAR       X       Y     FLUX     SHARP    ROUND Sky    SkyErr    Flux     FluxErr') do begin

     if strmatch(tmp,'Radius of aperture*') then begin
        tmp1=float(strmid(tmp,strpos(tmp,'=')+1))
        if napr eq 0 then aprsize=tmp1 else aprsize=[aprsize,tmp1]
        print,'Aperture ',napr,' has a size of',aprsize(napr)
        napr++
     endif
     readf,1,tmp
     nhlines++
  endwhile
  close,1



  nstars=(file_lines(subdir+'/'+filenames[0]+'.phot')-nhlines)*2/3
; if napr gt 6 then nstars=nstars/2	;there are two lines per star per polarization component if the number of apertures is greater than 6
  npairs=nstars/2

  numcol=8+napr*2
  tmp2=dblarr(numcol)

  starnum=intarr(nwp,nstars)
  pairnum=round(findgen(nstars)/2+0.5)
  xval=dblarr(nwp,nstars)
  yval=dblarr(nwp,nstars)
  sky=dblarr(nwp,nstars)
  skyerr=dblarr(nwp,nstars)
  flux=dblarr(nwp,napr,nstars)
  fluxerr=dblarr(nwp,napr,nstars)
  noise=dblarr(nwp,napr,nstars)
  sn=dblarr(nwp,napr,nstars)
  bestindex=intarr(nwp,nstars)
  bestaprid=intarr(nwp,nstars)
  psi=dblarr(nwp)
  z=dblarr(nwp)
  dz=dblarr(nwp)
  fluxo=dblarr(nwp)
  dfluxo=dblarr(nwp)
  fluxe=dblarr(nwp)
  dfluxe=dblarr(nwp)
  junk=''

  q2=dblarr(napr)
  dq2=dblarr(napr)
  u2=dblarr(napr)
  du2=dblarr(napr)
  p2=dblarr(napr)
  sigma2=dblarr(napr)
  dp2=dblarr(napr)
  theta2=dblarr(napr)
  dtheta2=dblarr(napr)

  z2=dblarr(nwp,napr)
  dz2=dblarr(nwp,napr)

  format='(12X,I5,2F8.2,F9.1,2F9.2,'+strtrim(napr*2+2,2)+'(F12.2))'


  for i=0,nwp-1 do begin
     openr,1,subdir+'/'+filenames[i]+'.phot'
     for k=0,nhlines-1 do readf,1,tmp
     if ~strcmp(tmp,'             STAR       X       Y     FLUX     SHARP    ROUND Sky    SkyErr    Flux     FluxErr') then begin
        print,'ERROR: '+filenames[i]+'.phot file does not have same number of header lines as '+filenames[0]+'.phot'
        print,'Check Shift between images'
        print,'...Quitting.'
        retall
     endif
	

	
	;loop through file reading in each line
     for j=0,nstars-1 do begin
        if ~(j MOD 2) then readf,1,junk ;skip every third line which lists the pair number
        readf,1,f=format,tmp2

        starnum(i,j)=tmp2(0)
        xval(i,j)=tmp2(1)
        yval(i,j)=tmp2(2)
        sky(i,j)=tmp2(6)
        skyerr(i,j)=tmp2(7)

        flux(i,*,j)=tmp2(8:7+napr)
        fluxerr(i,*,j)=tmp2(8+napr:7+2*napr)
        noise(i,*,j)=sqrt(flux[i,*,j]+!PI*aprsize(*)*aprsize(*)*(sky[i,j]+rdn^2))
        sn(i,*,j)=flux(i,*,j)/noise(i,*,j)
        maxindex=where(sn[i,*,j] gt [-9999,transpose(sn[i,0:napr-2,j])] and sn[i,*,j] gt [transpose(sn[i,1:napr-1,j]),-9999])
        bestindex(i,j)=maxindex(0)
     endfor
     close,1
     psi(i)=22.5*(i)*!PI/180
  endfor



; print,npairs
  openw,1,subdir+'/'+filenames[0]+'.pol'
  openw,4,subdir+'/'+filenames[0]+'all.pol'

  close,3
  openw,3,subdir+'/'+'calcpol.out'
  printf,3,'   ID 	X 	Y  	BestApp'

  reject_psigma=0
  reject_chisqrt=0
  nonconverges=0


  print,'Chisqrt limit =',Chisqrt_lim
  window,12,xsize=700,ysize=400
  plot,[0,100],[0,2*Chisqrt_lim],ytitle='Chi-square',xtitle='Polarization (%)',CharSize=1.5,/nodata
  oplot,[0,50,100],[Chisqrt_lim,Chisqrt_lim,Chisqrt_lim],linestyle=1
  XYOUTS,70,Chisqrt_lim, 'Chi-square limit', CHARSIZE=1.5, CHARTHICK=1

  for i=0,npairs-1 do begin


     printf,4,'PAIR',i+1
     printf,4,'             STAR       X       Y'
     printf,4,f='(12X,I5,2F8.2)',starnum[0,i*2],xval[0,i*2],yval[0,i*2]
     printf,4,f='(12X,I5,2F8.2)',starnum[0,i*2+1],xval[0,i*2+1],yval[0,i*2+1]


     oid=i*2
     eid=i*2+1
     sumo=0
     sume=0
     for k=0,napr-1 do begin
        for j=0,nwp-1 do begin
           bestaprid(j,i)=k     ;max([bestindex(j,oid),bestindex(j,eid)])
                                ;If bestaprid is -1, then set flux values to NAN
           if bestaprid(j,i) eq -1 then begin
              fluxo(j)=!Values.F_NAN
              dfluxo(j)=!Values.F_NAN
              fluxe(j)=!Values.F_NAN
              dfluxe(j)=!Values.F_NAN
           endif else begin
              fluxo(j)=flux(j,bestaprid(j,i),oid)
              dfluxo(j)=fluxerr(j,bestaprid(j,i),oid)
              fluxe(j)=flux(j,bestaprid(j,i),eid)
              dfluxe(j)=fluxerr(j,bestaprid(j,i),eid)
           endelse
	

        endfor
        sumo=total(fluxo)
        dsumo=sqrt(total(dfluxo*dfluxo))
        sume=total(fluxe)
        dsume=sqrt(total(dfluxe*dfluxe))
        ak=sume/sumo
        dak=sqrt(dsume*dsume+ak*ak*dsumo*dsumo)/abs(sumo)



        for j=0,nwp-1 do begin
           z(j)=(fluxe(j)-fluxo(j)*ak)/(fluxe(j)+fluxo(j)*ak)
           dz(j)=2.*sqrt((fluxo(j)*fluxo(j)*ak*ak*dfluxe(j)*dfluxe(j))+(fluxe(j)*fluxe(j)*ak*ak*dfluxo(j)*dfluxo(j))+(fluxe(j)*fluxe(j)*fluxo(j)*fluxo(j)*dak*dak))/((fluxe(j)+fluxo(j)*ak)* (fluxe(j)+fluxo(j)*ak))
        endfor

     
 
        q=total(z*round(cos(4.*psi)))/(nwp/2.)
        dq=(2./nwp)*sqrt(total(dz*dz*abs(round(cos(4*psi)))))
        u=total(z*round(sin(4.*psi)))/(nwp/2.)
        du=(2./nwp)*sqrt(total(dz*dz*abs(round(sin(4*psi)))))
        p = sqrt (q*q + u*u)
        dp=(1./p)*sqrt(q*q*dq*dq+u*u*du*du)
	
        sigma = sqrt ((total(z*z)/(nwp/2.) - q*q - u*u)/(nwp-2.))
        sigma2[k]= sigma
        sigma2_min=min(sigma2,Min_Subscript)

        
        theta=atan(u/q)

;        dtheta=0.5/(1.+(u/q)^2)*sqrt((du/q)^2+(-u*dq/q^2)^2)
;        dtheta=0.5/p^2*sqrt(q^2*du^2 + u^2*dq^2) ;(Naghizadeh-Khovei & Clarke 1993)
        dtheta=sigma/(2*p) ;assuming sigmau=sigmaq=sigmaP (Naghizadeh-Khovei & Clarke 1993)

        

        if q lt 0. then theta = theta + !PI
        if u lt 0. and q gt 0. then theta = theta + 2.*!PI
        theta = theta/2. 
        if theta ge !PI then theta = theta - !PI
        theta = !PI - theta + deltatheta *!PI/180. ;;>~2005

        if theta ge !PI then theta = theta - !PI
     
        q = p*cos(2.*theta)
        u = p*sin(2.*theta)


        theta=theta*180./!PI
        dtheta=dtheta*180./!PI

        z2[*,k]=z
        dz2[*,k]=dz
        
        q2[k]=q
        dq2[k]=dq
        u2[k]=u
        du2[k]=du
        p2[k]=p
        sigma2[k]=sigma
        dp2[k]=dp
        theta2[k]=theta
        dtheta2[k]=dtheta



;;           printf,1,' Best aperture size for each waveplate position  '
        printf,4,f='(2x,A,'+strtrim(nwp,2)+'I5)','Waveplate Position:',findgen(nwp)+1
        printf,4,f='(2x,A,'+strtrim(nwp,2)+'I5)','Aperture Size:',round(aprsize(bestaprid(*,i)))
        printf,4,f='(2x,A,'+strtrim(nwp,2)+'F12.5)','z values:          ',z
        printf,4,f='(2x,A,'+strtrim(nwp,2)+'F12.5)','dz values:          ',dz
        printf,4,'     q       dq        u       du        p     sigma      dp     theta   dtheta'
        printf,4,f='(7(F9.5," "),2(F8.3," "))',q,dq,u,du,p,sigma,dp,theta,dtheta
              

     endfor                     ;naper; loop over k     

     printf,4,f='(2x,A,'+strtrim(nwp,2)+'I5)','Best Aperture (aperture with lowest psigma):',Min_Subscript+1
     printf,4,''

 

     for k=0,napr-1 do begin
        if (k eq  Min_Subscript) then begin
 
                 printf,3,f='(I5,2(F8.2),I5)',starnum[0,i*2],xval[0,i*2],yval[0,i*2],Min_Subscript+1
                 printf,3,f='(I5,2(F8.2),I5)',starnum[0,i*2+1],xval[0,i*2+1],yval[0,i*2+1],Min_Subscript+1
          
           Chisqrt=graf(z2[*,k],dz2[*,k],pair=i+1)     ;Returns reduced Chisqrt to discriminate bad zfits
           if FINITE(Chisqrt) eq 0 then nonconverges++ ;counting the non converge fits


           plots,p2[k]*100.,chisqrt,PSYM=1 ;plotting chisqrt vs pol

           if (Chisqrt gt Chisqrt_lim) then begin   
              print,'Pair rejected by Chi-square: ',i+1,', Chisqrt = ',Chisqrt
              reject_chisqrt++
           endif else begin   


              if (p2[k]/sigma2[k] ge p_over_sigma*1.)  and (p2[k] le 1.) then begin


                 printf,1,'PAIR',i+1
                 printf,1,'             STAR       X       Y'
                 printf,1,f='(12X,I5,2F8.2)',starnum[0,i*2],xval[0,i*2],yval[0,i*2]
                 printf,1,f='(12X,I5,2F8.2)',starnum[0,i*2+1],xval[0,i*2+1],yval[0,i*2+1]
                 printf,1,f='(2x,A,'+strtrim(nwp,2)+'I5)','Waveplate Position:',findgen(nwp)+1
                 printf,1,f='(2x,A,'+strtrim(nwp,2)+'I5)','Best Aperture Size:', round(dblarr(nwp)+k+1)
                 printf,1,f='(2x,A,'+strtrim(nwp,2)+'F12.5)','z values:          ',z2[*,k]
                 printf,1,f='(2x,A,'+strtrim(nwp,2)+'F12.5)','dz values:          ',dz2[*,k]
                 printf,1,'     q       dq        u       du        p     sigma      dp     theta   dtheta'
                 printf,1,f='(7(F9.5," "),2(F8.3," "))',q2[k],dq2[k],u2[k],du2[k],p2[k],sigma2[k],dp2[k],theta2[k],dtheta2[k]
                 printf,1,''



              endif else begin  ;p/s
                 if FINITE(p2[k]/sigma2[k]) eq 1 then begin 
                    print,'Pair rejected by p/sigma (or spurous pol greater that 100%)',i+1,', P%=',p2[k]*100.,'+/-',sigma2[k]*100 
                    reject_psigma++
                 endif
              endelse
           endelse              ;chisqrt
        endif
     endfor 

  endfor                        ;pairs; loop over i   
  close,1
  close,3
  close,4

  print,''
  print,'File with all the stars before p/sigma rejection: '+subdir+'/'+filenames[0]+'all.pol'
  print,'File with the stars after p/sigma rejection: '+subdir+'/'+filenames[0]+'.pol'
  print,''

  print,'No. of Pairs rejected by p/sigma creteria (or spurous pol greater that 100%): ',reject_psigma
  print,'No. of Pairs rejected by failed to converge modulation fit: ',nonconverges
  print,'No. of Pairs rejected by Chi-square creteria: ',reject_chisqrt
  print,'Chisqrt limit > ',Chisqrt_lim

if KEYWORD_SET(chisqrt_int) eq 1 then begin 
  Chisqrt_lim = ' '             ;chisqrt limit to reject bad modulation fit
  print,'Would you like to change the Chisqrt limit value?'
  read,' enter the value, else press enter to continue : ', Chisqrt_lim

  IF (Chisqrt_lim ne '') THEN GOTO, JUMP1
endif

end
