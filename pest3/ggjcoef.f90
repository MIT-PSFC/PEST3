  SUBROUTINE pstggjcoef(isurf,pe,pf,ph,pdi,pdr,pak,pgee,pem,&
                        feq,zor2,zvp,pq0s,pls)
!....................................................................
!
! compute the glasser et al. coef. at surface isurf
!
! aplet  10.1_r8 .97
!
  USE pstcom
  USE mtrik1
  USE newmet
  USE l22com
  USE temps
  USE r33com
  USE comggf
  IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER ISURF
      real*8 PE
      real*8 PF
      real*8 PH
      real*8 PDI
      real*8 PDR
      real*8 PAK
      real*8 PGEE
      real*8 PEM
      real*8 FEQ
      real*8 PQ0S
      real*8 PLS
      INTEGER NR
      real*8 G2
      real*8 PP
      real*8 QP
      real*8 GG
      real*8 ZNQP2
      real*8 ZVD
      real*8 ZVP
      real*8 ZVPP
      real*8 ZOR2
      real*8 ZB2
      real*8 ZBOG2
      real*8 ZB2OG2
      real*8 ZB2G2
      real*8 ZOB2G2
      real*8 ZOG2
      real*8 ZOB2
      real*8 ZOB3
      real*8 ZG2OB2
      real*8 ZG2OB3
      real*8 ZB
      real*8 ZALPH
      real*8 ZLAMBD
      real*8 ZRADS
      real*8 Z
      real*8 DX
      real*8 dMdB

      INTEGER :: t,i
      REAL*8, dimension(nths) :: zbsq
 
      nr = isurf
!
! new coefficients by aplet  7.1_r8 .93
!
!
      g2 = r*r * ga(nr)**2
      pp = ppa(nr)
      qp = qpa(nr)
      gg = ga(nr)*r
      znqp2 = (n * qpa(nr))**2
 
      zvp =  0._r8 
      zvpp =  0._r8 
      zb2 =  0._r8 
      zb2og2 =  0._r8 
      zbog2 =  0._r8 
      zb2g2 =  0._r8 
      zob2g2 =  0._r8 
      zog2 =  0._r8 
      zob2 =  0._r8  
      zob3 =  0._r8  
      zg2ob2 =  0._r8 
      zg2ob3 =  0._r8 
      zb =  0._r8 
      zalph =  0._r8 
      zlambd =  0._r8 
      zrads =  0._r8 

      zbsq(1:mth) = ( grpssq(1:mth,nr)+g2 )/ xsq(1:mth,nr)
 
      dth = twopi / mth

      zvp    = sum( xjacob(1:mth,nr)                                )*dth
      zvd=0

      do i=1,nosurf-1
        dx=psia(i+1)-psia(i)
        zvd = zvd+ (sum(xjacob(1:mth,i))*dth &
                  + sum(xjacob(1:mth,i+1))*dth)/2  *dx
      enddo

      ! xjprym = d xjacob / d Psi; Psi = psi/twopi

zvpp   = sum( xjprym(1:mth,nr)                                )*dth/twopi
zb     = sum( xjacob(1:mth,nr)*sqrt(zbsq(1:mth))              )*dth
zb2    = sum( xjacob(1:mth,nr)*zbsq(1:mth)                    )*dth
zbog2  = sum( xjacob(1:mth,nr)*sqrt(zbsq(1:mth))/grpssq(1:mth,nr) )*dth
zb2og2 = sum( xjacob(1:mth,nr)*zbsq(1:mth)/grpssq(1:mth,nr)   )*dth
zob2g2 = sum( xjacob(1:mth,nr)/(zbsq(1:mth)*grpssq(1:mth,nr)) )*dth
zog2   = sum( xjacob(1:mth,nr)/grpssq(1:mth,nr)               )*dth
zob2   = sum( xjacob(1:mth,nr)/zbsq(1:mth)                    )*dth
zob3   = sum( xjacob(1:mth,nr)/zbsq(1:mth)**1.5               )*dth
zor2   = sum( xjacob(1:mth,nr)/xsq(1:mth,nr)                  )*dth
zg2ob2 = sum( xjacob(1:mth,nr)*grpssq(1:mth,nr)/zbsq(1:mth)   )*dth
zg2ob3 = sum( xjacob(1:mth,nr)*grpssq(1:mth,nr)/zbsq(1:mth)**1.5 )*dth

      zalph  = twopi * n * zvp
      zlambd = - qp/( twopi*zvp**3 )
 
      pe   = - pp*zb2og2*(zvpp - gg*qp*zvp/zb2)/(twopi*qp**2)
      pf   = pp**2 *(zb2og2*zob2 + g2*(zb2og2*zob2g2 - zog2**2)) &
             /(twopi*qp)**2
      ph   = gg*pp*(zog2 - zvp*zb2og2/zb2)/(twopi*qp)
      pdi = pe + pf + ph -  0.25_r8 
      pdr = pe + pf + ph*ph
 
 
      pem  = zb2og2*( zg2ob2 + g2*(zob2 - zvp**2/zb2) )/zvp**2
      dMdB = 2.*zbog2*( zg2ob2 + g2*(zob2 - zvp**2/zb2) )/zvp**2  &
              + zb2og2*(-2.*zg2ob3+2.*g2*(-zob3+zvp**2*zb/zb2**2))/zvp**2
!
! pq0s and pls are the typical growth rate factors and layer width divided by
! eta^1/3
! 
!     write(*,*) 'zb',zb
!     write(*,*) 'zb2',zb2
!     write(*,*) 'zob2',zob2
!     write(*,*) 'zob3',zob3
!     write(*,*) 'zbog2',zbog2
!     write(*,*) 'zb2og2',zb2og2
!     write(*,*) 'zg2ob2',zg2ob2
!     write(*,*) 'zg2ob3',zg2ob3
!     write(*,*) nr
!     write(*,*) 'zb2og2',zb2og2
!     write(*,*) 'zvp',zvp
!     write(*,*) 'twopi',twopi
!     write(*,*) 'znqp2',znqp2
!     write(*,*) 'qp',qp

      zor2=zor2/zvp

!     write(*,*) '1/<r2>',zor2
!     write(*,*) 'f', sqrt(g2)
      feq=sqrt(g2)

      z   =  zb2/zb2og2
      pq0s = abs( (twopi**2 *znqp2 *z)/(zvp**2 *pem )  &
          )**( 1._r8 / 3._r8 )
      pls = abs( (zvp**2 *pem *z**2)/(twopi**2 *znqp2)  &
          )**( 1._r8 / 6._r8 )

!     write(*,*) 'pq0s',pq0s
!     write(*,*) 'pls',pls
!     write(*,*) 't_a',(zvp/twopi)*(pem/znqp2)**0.5_r8
!     write(*,*) 't_r',1._r8/z
!     write(*,*) 'pq0s*t_a=pls',pq0s*(zvp/twopi)*(pem/znqp2)**0.5_r8
!     write(*,*) 'elecd*s',1._r8/(z*(zvp/twopi)*(pem/znqp2)**0.5_r8)
!     write(*,*) 'dMdB',dMdB
 
      if( abs(pp) <  1.E-6_r8  ) then
         pak =  9.E99_r8 
      else
         pak  = (twopi**2 *zvp**2 *zlambd)**2 *z/(pem* pp**2)
      end if
      if( abs( pa(nr) ) <  1.E-6_r8  ) then
         pgee =  9.E99_r8 
      else
         pgee = zb2/(gamma*pa(nr)*pem*zvp)
      end if    

!     write(*,*) 'pe',pe
!     write(*,*) 'pf',pf
!     write(*,*) 'pgee',pgee
!     write(*,*) 'ph',ph
!     write(*,*) 'pak',pak
!     write(*,*) 'pem',pem

!     Write out the inner layer code input file.  

!          write(outmod, 1480) i, xinf(i), zinf(i), xwal(i), zwal(i)
!
!
!1470 format(//' plasma and wall boundaries:' &
!/ 1x,' i ',1x,5x,'xinf',  11x,'zinf',11x,'xwal',11x,'zwal')
!
!1480 format(1x,i3,4(1x,e14.7)) 

      return
 
      end


