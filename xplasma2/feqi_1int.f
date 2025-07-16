      REAL*8 function feqi_1int(iorder,zrho,zspl,nrho,j,zdelx)
c
c  compute integral from zrho(j) to zrho(j)+zdelx
c  zrho(j)+zdelx .le. zrho(j+1) expected.
c
      use xdistrib_mod
c
      IMPLICIT NONE
c
c  input:
      integer iorder                    ! order (1:Hermite, 2:Spline)
      integer nrho                      ! no. of pts.
      REAL*8 zrho(nrho)                 ! indep coord.
      REAL*8 zspl(ncdim(iorder),nrho) ! spline of function
c
      integer j                         ! index to zone
      REAL*8 zdelx                      ! displacement w/in the zone
c
c-------------------------
c
      REAL*8 zh,zhi,zh2,zh3
c
      REAL*8 f0,f1,fx0,fx1,fxx0,fxx1
      REAL*8 p,p2
      REAL*8 zterm1,zterm2,zterm3,zterm4,zsum
c
      real*8, parameter :: cone=1,ctwo=2,cthree=3,cfour=4,csix=6
      real*8, parameter :: chalf=1/ctwo, cthird=1/cthree, c4th=1/cfour
c
c-------------------------
c
      f0=zspl(1,j)
      if((iorder.ge.0).and.(iorder.le.2)) then
         zh=zrho(j+1)-zrho(j)
         zhi=cone/zh
         p=zhi*zdelx
         p2=p*p
         f1=zspl(1,j+1)
      endif
c
      if(iorder.eq.3) then
c
         zsum=zdelx*(f0+ chalf*zdelx*(zspl(2,j) +
     >      cthird*zdelx*(zspl(3,j) + c4th*zdelx*zspl(4,j))))
c
c  explicit spline formulation
c
      else if(iorder.eq.2) then
c
c  compact spline formulation
c
         zh3=zh*zh*zh
c
         fxx0=zspl(1,j)
         fxx1=zspl(1,j+1)
c
         zterm1=zh*f0*(p-p2/ctwo)
         zterm2=zh*f1*(p2/ctwo)
         zterm3=zh3*fxx0*p2*(-cone+p-p2/cfour)/csix
         zterm4=zh3*fxx1*p2*(p2/cfour-cone/ctwo)/csix
         zsum=zterm1+zterm2+zterm3+zterm4
c
      else if(iorder.eq.1) then
c
c  Hermite formulation
c
         zh2=zh*zh
c
         fx0=zspl(1,j)
         fx1=zspl(1,j+1)
c
         zterm1=zh*f0*(p2*(p2/ctwo-p)+p)
         zterm2=zh*f1*(p2*(p-p2/ctwo))
         zterm3=zh2*fx0*p2*(p2/cfour-ctwo*p/cthree+chalf)
         zterm4=zh2*fx1*p2*(p2/cfour-p/cthree)
         zsum=zterm1+zterm2+zterm3+zterm4
c
      else if(iorder.eq.0) then
c
         zsum=zh*p*(f0 + p*(f1-f0)/ctwo)
c
      else if(iorder.eq.-1) then
c
         zsum=zdelx*f0
c
      else
         write(6,'(a)') ' ?eqi_1int:  iorder out of range!'  ! RGA:dep remove eq_errmsg()
      endif
c
      feqi_1int=zsum
c
      return
      end
