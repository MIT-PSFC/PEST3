#ifdef __DEBUG

subroutine scrdesdbg(zlabel,rmb,zmb,nmb,nmbout,rmbout,zmbout)

  ! descur/desreg debug plotting tool
  ! plot vs. theta of: arclength, dist from (R0,Z0), arctan((Z-Z0)/(R-R0))
  ! for both moments and for the original data in which "theta" is defined
  ! as twopi*(i-1)/(nmb-1)

  implicit NONE
  character*(*) zlabel
  integer nmb
  real*8 rmb(nmb),zmb(nmb)
  integer nmbout
  real*8 rmbout(0:nmbout,2),zmbout(0:nmbout,2)

  !-----------
  !  automatic...
  real*8 :: thmb(nmb),arclen(nmb),arcleno(nmb)
  real*8 :: dist(nmb),disto(nmb),angtan(nmb),angtano(nmb)
  real*8 :: zsin(nmbout),zcos(nmbout)
  character*60 mnlab(3)
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  real*8 :: zr0,zz0,zdl,zcorr
  real*8 :: zrcur,zzcur,zrpre,zzpre
  integer :: i,j
  !-----------

  mnlab=' '
  mnlab(1) = zlabel//' debug plot '

  ! from input data construct arclength and [R0,Z0]

  zr0=0; zz0=0
  do i=1,nmb
     thmb(i)=(i-1)*c2pi/(nmb-1)
     if(i.eq.1) then
        arclen(i)=0
     else
        zdl=sqrt((rmb(i)-rmb(i-1))**2+(zmb(i)-zmb(i-1))**2)
        arclen(i)=arclen(i-1)+zdl
        zr0=zr0+zdl*(rmb(i)+rmb(i-1))/2
        zz0=zz0+zdl*(zmb(i)+zmb(i-1))/2
     endif
  enddo
  zr0=zr0/arclen(nmb)
  zz0=zz0/arclen(nmb)

  ! dist. vs. theta & arctan vs. theta for input data

  zcorr=0
  do i=1,nmb
     dist(i)=sqrt((rmb(i)-zr0)**2+(zmb(i)-zz0)**2)
     angtan(i)=atan2(zmb(i)-zz0,rmb(i)-zr0) + zcorr
     if(i.eq.1) then
        if(angtan(1).lt.-c2pi/2) angtan(i)=angtan(i)+c2pi
        if(angtan(i).gt.c2pi/2) angtan(i)=angtan(i)-c2pi
     else
        if((angtan(i)-angtan(i-1)).gt.c2pi/2) then
           angtan(i)=angtan(i)-c2pi
           zcorr=zcorr-c2pi
        else if((angtan(i)-angtan(i-1)).lt.-c2pi/2) then
           angtan(i)=angtan(i)+c2pi
           zcorr=zcorr+c2pi
        endif
     endif
  enddo

  ! same stuff for the output data

  zcorr=0
  do i=1,nmb
     if(i.gt.1) then
        zrpre=zrcur
        zzpre=zzcur
     endif
     call r8sincos(thmb(i),nmbout,zsin,zcos)
     zrcur=rmbout(0,1)
     zzcur=zmbout(0,1)
     do j=1,nmbout
        zrcur=zrcur+rmbout(j,1)*zcos(j)+rmbout(j,2)*zsin(j)
        zzcur=zzcur+zmbout(j,1)*zcos(j)+zmbout(j,2)*zsin(j)
     enddo
     disto(i)=sqrt((zrcur-zr0)**2+(zzcur-zz0)**2)
     angtano(i)=atan2(zzcur-zz0,zrcur-zr0) + zcorr
     if(i.eq.1) then
        if(angtano(1).lt.-c2pi/2) angtano(i)=angtano(i)+c2pi
        if(angtano(i).gt.c2pi/2) angtano(i)=angtano(i)-c2pi
        arcleno(i)=0
     else
        if((angtano(i)-angtano(i-1)).gt.c2pi/2) then
           angtano(i)=angtano(i)-c2pi
           zcorr=zcorr-c2pi
        else if((angtano(i)-angtano(i-1)).lt.-c2pi/2) then
           angtano(i)=angtano(i)+c2pi
           zcorr=zcorr+c2pi
        endif
        zdl=sqrt((zrcur-zrpre)**2+(zzcur-zzpre)**2)
        arcleno(i)=arcleno(i-1)+zdl
     endif
  enddo

  mnlab(3)='arclength vs. theta'
  call r8_grafx2(thmb,arclen,arcleno,nmb,' ','L',mnlab(1),mnlab(2),mnlab(3))

  mnlab(3)='distance from [R0,Z0] vs. theta'
  call r8_grafx2(thmb,dist,disto,nmb,' ','L',mnlab(1),mnlab(2),mnlab(3))

  mnlab(3)='arctan to [R0,Z0] vs. theta'
  call r8_grafx2(thmb,angtan,angtano,nmb,' ','L',mnlab(1),mnlab(2),mnlab(3))

end subroutine scrdesdbg

#else

subroutine scrdesdbg

   write(6,*) ' %NOTE: dummy scrdesdbg call: compiled w/o __DEBUG '
   return

end subroutine scrdesdbg

#endif
