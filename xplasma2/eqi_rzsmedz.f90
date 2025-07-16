subroutine eqi_rzsmedZ(R,inumR,Z,inumZ,fdata,zsm)

  !  smooth in Z direction, in edge vicinity (see eqi_rzsmedg.for).

  use xplasma_definitions
  use eqi_rzbox_module

  implicit NONE

  integer inumR,inumZ
  real*8 R(inumR),Z(inumZ)
  real*8 fdata(inumR,inumZ)
  real*8 zsm

  !-----------------------------

  real*8 eps(inumZ),del(inumZ),f(inumZ),fsm(inumZ),zwkR(inumZ)
  real*8 dist(inumZ)

  real*8 zdista,zdum,zdum2

  integer iR,iZ,ict,idum1,idum2,iertmp

  real*8, parameter :: HALF = 0.5d0
  real*8, parameter :: C3OV2= 1.5d0
  real*8, parameter :: zbc = HALF

  !-----------------------------

  eps=1.0d30
  del=0.0d0

  do iR=1,inumR
     zwkR=R(iR)
     call xplasma_bdfind(sp,inumZ,ZwkR,Z,iertmp, &
          maptype=3,dist=dist)

     ict=0
     do iZ=1,inumZ
        zdista=abs(dist(iZ))
        if(zdista.le.HALF*zsm) then
           ict=ict+1
           del(iZ)=zsm
        else if(zdista.lt.C3OV2*zsm) then
           ict=ict+1
           del(iZ)=C3OV2*zsm-zdista
        endif
     enddo

     if(ict.gt.0) then

  !  something to smooth

        f=fdata(iR,1:inumZ)
        call r8filtr6(Z,f,fsm,inumZ,eps,inumZ,eps,0,del, &
             0,zdum,0,zdum,zbc,zbc,zdum2,idum1,idum2)
        fdata(iR,1:inumZ)=fsm

     endif
  enddo

end subroutine eqi_rzsmedZ
