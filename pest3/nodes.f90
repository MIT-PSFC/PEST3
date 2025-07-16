!
!......................................................................
      SUBROUTINE pstnodes
!......................................................................
!       version of nodes to permit arbitrary distribution of fin. elements
!       the plasma is divided into regions by the singular surfaces.
!       maximum allowed is 5 regions ie. 4 sing. surfaces.
!       each region is divided into zones where the nodeds are either
!       distributed uniformly, or packed with a power law around the
!       singular surfaces.
!
!       modified version of nodes to allow packing of mesh at several
!        resonant surfaces. mesh distribution is specified by:
!       alfmsh(5)  exponent of mesh scaling around magnetic axis and,alfmsh(i+1)
!                  rational surfaces.
!       widmsh(5)  width of scaled zone.
!       distribution function of mesh nodes is
!  -/+zg(i)*|(z-zts(i))/widmsh(i)|**alfmsh(i)+zd(i)    in scaled zone.
!       zb*z + zc(i)       in linear zone,
!       for z = inode-1/m,
!       and i=1,nsing+ 1._r8   zts, zg, zd, zb and zc are computed
!       to ensure continuity of distribution function and its
!       derivative. a. pletzer, 2 july  90._r8 
!
! mofifications:
!
!         uniform distribution of z within each section. aplet  20.4_r8 .93
!         introduce distortion function to increase node density
!         near plasma edge.                              aplet  6.5_r8 .93
! =======
! packmth
!    0    optional edge packing using standard zone DPB 01/10
!    1    optional edge packing with a different power law than standard
!           zone.  The point is to resolve the peeling mode very near the
!           separatrix.  The gradient in the zone is specified as
!
!           dzx/dz=a*(c0-zx)^alfa 
!
!           where a and c0 are found and alfa and w are set the same as
!           standard zones.  The continuity of the gradient to the linear
!           zone gives
!
!           a=b/(c0+wx-1)^alfa
!
!           Where wx is the width of the edge zone in zx space.  The gradient
!           is specified to be a set fraction f of the background b at zx=1.  
!           This gives c0 as 
!
!           c0=(fb/a)^(1/alfa)+1=(1+f^(1/alfa)(wx-1))/(1-f^(1/alfa))
!
!           The function integrates to 
!
!           zx=c0-((1-alfa)(-a*z - c1))^(1/(1-alfa)) 
!
!           Thus alfa is constrained to alfa>1. The boundary on zx at z=1, 
!           and the continuity at zx=1 give an analytic expression for c1
!
!           c1=-(d-1)^(1-alfa)/(1-alfa)-a
!
!           and an equation whos root is found numerically using zridderx 
!           for b, respectively. This produces a mesh that is packed near
!           the edge with many points as the gradient is small and the packing
!           extends inside of the edge.  This is inspired by mesh requirements
!           with GATO and step lengths with DCON in resolving the (stable) 
!           peeling modes, but is not likely necessary for equilibria without
!           significant edge current.  DPB 2/10
!
!    2    grid packing is changed to a smoothly transitioning packing 
!           function 
!
!           dzx/dz=((z-z1)/(z2-z1))^alfa*((z2-z)/(z2-z1))^beta
!
!           where z1 is the lower bound z of a region between surfaces (or
!           axis or edge) and z2 is the upper bound z of the next surface
!           etc. alfa and beta are the given packing exponents at each 
!           surface respectively.  Here, since we want the grid distribution
!           proportional to 1/mu, and its derivative proportional to 
!           1/mu-1, the alfa and beta are set to alfmsh(i)-1 and
!           alfmsh(i+1)-1 for surfaces i and i+1 respectively.  The -1
!           is thus absorbed in the actual code below, compared to these 
!           notes. 
!
!           Defining t=(z-z1)/(z2-z1) we have
!
!           dzx/dz=t^alfa*(1-t)^beta 
!
!           which integrates to the hypergeometric function 2F1 for zx
!
!           zx=C*(t^(1+alfa)/(1+alfa))*2F1(1+alfa,-beta,2+alfa,t)
!
!           C here is the standard normalization constant 
!
!           C=Gamma(2+alfa+beta)/(Gamma(1+alfa)*Gamma(1+beta))
!
!           A numerical recipes based algorithm is implemented to calculate
!           this grid distribution.  The intention is to have the appropriate
!           power in the grid density in the limits at t=0 and t=1 and 
!           avoid having a background grid density, which will allow the 
!           grid to become more sparse in the inter-surface zones where 
!           high grid density is not needed.  This will allow us to access
!           more highly packed near surface grids, and resolve cases 
!           otherwise unresolvable due to memory and machine limitations.
!......................
      USE pstcom

      USE temps

      USE comggf

      USE nrtype

      USE hypgeo_info

      USE nr, only : beta,hypgeo

      IMPLICIT NONE
      complex(DPC) :: ahg,bhg,chg
      complex(DPC), dimension(:), allocatable :: xhg,dz,z1hg,z2hg,z3hg,&
                                                 z4hg,z5hg,z6hg
      real(DP) :: normfac

      integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)
      INTEGER 		:: N2
      INTEGER 		:: MTOT
      INTEGER 		:: JTEST
      INTEGER 		:: J
      INTEGER 		:: JMAX
      INTEGER 		:: J2
      INTEGER 		:: MMINUS
      INTEGER 		:: I
      INTEGER 		:: JFLAG2
      INTEGER 		:: JFMAX
      REAL*8 		:: SUM
      REAL*8 		:: B
      REAL*8 		:: Z
      INTEGER 		:: JFLAG
      REAL*8 		:: Z1
      REAL*8 		:: Z2
      INTEGER 		:: I3
      INTEGER 		:: I0
      INTEGER 		:: I1
      INTEGER 		:: I2
      REAL*8 		:: ZM
      INTEGER 		:: JS
      REAL*8 		:: ZDXS
      INTEGER 		:: IN
      INTEGER 		:: IN2
      REAL*8 		:: ZX
      REAL*8 		:: DELPSL
      REAL*8 		:: DELPSR
      REAL*8 		:: DELPSI
      REAL*8 		:: DELPSA
      INTEGER 		:: KN
      INTEGER 		:: KM
      REAL*8 		:: DPS
      INTEGER 		:: IK
      REAL*8 		:: ZL
      REAL*8 		:: ZR
 
      LOGICAL	 	 	 :: nodchk
      REAL*8, DIMENSION(nfe) 	 :: width
      REAL*8, DIMENSION(nsg1) 	 :: zc
      REAL*8, DIMENSION(nsg1) 	 :: zg
      REAL*8, DIMENSION(nsg1) 	 :: zd
      REAL*8, DIMENSION(nsg1) 	 :: zts
      REAL*8, DIMENSION(nsg1) 	 :: zsing
      REAL*8	 	 	 :: cf
      REAL*8	 	 	 :: cs
      REAL*8	 	 	 :: arad
      REAL*8, DIMENSION(nsg1,nsg1) 	 :: zonew
      REAL*8, DIMENSION(nsg) 	 :: xrat
      INTEGER, DIMENSION(55) 	 :: iord
      INTEGER, DIMENSION(5) 	 :: lzone
      INTEGER, DIMENSION(2) 	 :: iop
      INTEGER, DIMENSION(3,nsg1) 	 :: mzoned
      INTEGER, DIMENSION(2) 	 :: minner
      INTEGER	 	 	 :: pkptsum

      INTEGER 		:: index_ratsurf
      REAL*8 		:: zdiav

      LOGICAL 		:: iok(1) ! if iok(i) TRUE then SKIP element i of vector
      REAL*8 		:: x(1)     ! input vector {x(i)}
      REAL*8 		:: ansr(1)  ! output vector {f(x(i))}
      INTEGER 		:: ninput,noutput
      REAL*8 		:: zinput(1,4),zoutput(1,1)

      REAL*8 		:: r1,b1,b2,res,oldres,bs,we,af,tol,wx,c0,c1
      REAL*8 		:: jfac,sumfac
      INTEGER 		:: found,k,ifail,pkmth,dlaymth,shft

      COMMON /c2f90nodchk/ nodchk 
      COMMON /c2f90zonew/ zonew 
      COMMON /c2f90mzoned/ mzoned 
      COMMON /c2f90minner/ minner 
      COMMON /c2f90arad/ arad 
      COMMON /c2f90xrat/ xrat 
      EXTERNAL resid_b

      dlaymth=1

      if( imode == 1 ) then
        write(itty,8756) nosurf,dlayb, dlay
        write(outmod,8756) nosurf,dlayb, dlay
 8756   format(/' equilibrium grid size         nosurf = ',i7/  &
        ' rel. width of big element      dlayb = ',f13.6,/  &
        ' rel. width for frobenius coeff. dlay = ',f13.6/)
      end if
!
      write(itty,8755) m
      write(outmod,8755 ) m
 8755 format(/' number of finite elements          m = ',i7/)
!
      msing(1) = 1
      nusurf   = nosurf
      mtot = m
!
      psidif = psia(nosurf) - psia(1)

!     take out the local symmetric packing points around each rational surface
      mminus=m
      do i=1,nosing
        mminus = mminus - 2*packpts(i)
      enddo
!
! test for presence of singular surfaces...
!
      if( nosing ==  0 ) then
         !
         ! linear mesh...
         !
         psinod = (/ (psidif*(i-1)/m, i=1,mtot+1) /) + psia(1)
         !
         go to 101
      end if
!
! there are singular surface(s)...

! this sets widmsh=0 for the first zone outside the packed grid 
      widmsh(nsing1+2) =  0.0_r8 

! compute mesh exponents
      
! the alfmsh of the axis and edge with the optional edge packing is handled 
! on start up.  With alfmsh=1 in any zone there is no packing, this is the
! default for the edge. To avoid packing the surfaces use lpack.  DPB

!     lpack=.false.
      if(lpack) then
         if (nsing1>=2) then
            do i = 2, nsing1 
               if(abs(xmu(i-1)-0.5_r8)>1.e-6_r8) alfmsh(i) = 1.0_r8/xmu(i-1)
            enddo
         endif
      endif

!
!...............................................
!
       jflag2 = 0
       jfmax = 10
!
 100   continue

       if (packmth/=2) then
!
!       check alfmsh and widmsh...
!
         sum = - widmsh(1) - widmsh(nsing1+1)
       
         do 80 i=1,nsing1+1
         alfmsh(i) = max(1.0_r8, abs( alfmsh(i) ) )
         widmsh(i) = abs( widmsh(i) )
         sum = sum +  2._r8 *widmsh(i)
 80       continue
!
         if( sum >=    1._r8  ) then
         write(itty,*)  ' *** widmsh too large'
         write(outmod,*)' *** widmsh too large'

         do 81 i=1,nsing1+1
         widmsh(i) = widmsh(i)/(sum+ 0.1_r8 )
 81      continue
         write(itty,*)'corrected values are:'
         write(outmod,*)'corrected values are:'
         write(itty,*)'widmsh = ',(widmsh(i),i=1,nsing1+1) 
         write(outmod,*)'widmsh = ',(widmsh(i),i=1,nsing1+1) 
         end if
       end if
!
!       position of rational surfaces...
!
       do 82 i=1,nsing1+1
       zsing(i) = (psisin(i) - psia(1))/psidif
 82    continue
!
       if (packmth/=2) then
!      slope of mesh distribution function in linear zones...
!      
!      insert new code for b here with option on packmth
!
!      first add up the z widths of the zones leading up to the edge.
!
         b =  1._r8  - (alfmsh(1)- 1._r8 )*widmsh(1)/alfmsh(1)
         do i=2,nsing1
            b = b -  2._r8 *(alfmsh(i)- 1._r8 )*widmsh(i)/alfmsh(i)
         end do
!      
         select case(packmth)
           case(0)
!            then, if the last zone is standard
!            add the last regular (edge) zone into the width
             b =  b  - (alfmsh(nsing1+1)- 1._r8 )* &
                           widmsh(nsing1+1)/alfmsh(nsing1+1)
!            the inverse is b.
             b =  1._r8 /b
           case(1)
!            find the b numerically; at this point, b variable is 
!            the width of all zones minus the edge zone (not inverted)
             cs = b-widmsh(nsing1+1)
!            the fraction of b in the edge packing at the boundary 
             cf=0.01_r8
!            this is somewhat simple.  Here I brackett the first root for b by
!            stepping through b space from 0 to 3.  If b is bigger than 3 it
!            will fail. If so, simply adjust.
             res=0
             we=widmsh(nsing1+1)
             af=alfmsh(nsing1+1)
             ninput=4
             zinput(1,1)=cf
             zinput(1,2)=we
             zinput(1,3)=cs
             zinput(1,4)=af
             noutput=1
             found=0
             do k=1,3*m
               oldres=res
               bs=real(k)/real(m)
               call resid_b(1,.true.,bs,res, &
                            1,zinput,ninput,zoutput,1)
               if (k>1 .and. found==0 .and. res*oldres<0) then
                 found=1
                 r1=res
                 b2=bs
                 b1=real(k-1)/real(m)
               endif
             enddo
!            find the root bracketted by b1 b2 put into b
             tol = 1.0e-12_r8
             ifail=0
             call zridderx(1,.false.,b1, b2, tol, tol*r1, &
                  resid_b, b, ifail, 1, zinput, ninput, zoutput, 1)
  
             if(ifail.ne.0) then
               stop 'ERROR in nodes: Error in root finding for b' 
             end if
           
             wx=1-b*cs
             c0=(1+cf**(1/af)*(wx-1))/(1-cf**(1/af));
             c1=-(c0-1)**(1-af)/(1-af)-b/(c0+wx-1)**af
  
           case default
             print *, 'packmth ',packmth,' not handled'
             stop 'ERROR in nodes: packmth not handled'
         end select

!
!         first zone near magnetic axis...
!
         zd(1) =  0._r8 
         zg(1) = b*widmsh(1)/alfmsh(1)
!
!         first linear mesh zone...
!
         zc(1) = - b*widmsh(1)*(alfmsh(1)-1)/alfmsh(1)
!
!         iterate over regions...
!
         zts(1) =  0._r8 
         do i=2,nsing1
            zc(i) = zc(i-1) -  2._r8 *b*(alfmsh(i)-1)*widmsh(i)/alfmsh(i)
            zg(i) = b*widmsh(i)/alfmsh(i)
            zd(i) = zsing(i)
          !
          !aplet  22.4_r8 .98 zts(i) = ( zd(i) - (zc(i)+zc(i-1))/ 2._r8  )/b
            zts(i) = zts(i-1) - zsing(i-1)/b + zsing(i)/b &
                 + widmsh(i-1)*(alfmsh(i-1)- 1._r8 )/alfmsh(i-1) &
                 + widmsh(i  )*(alfmsh(i  )- 1._r8 )/alfmsh(i  )
         end do

!        The last zones at the edge.

         select case(packmth)
           case(0)
!            zts(nsing1+1) = zts(nsing1)-zsing(nsing1)/b+zsing(nsing1+1)/b &
!              + widmsh(nsing1)*(alfmsh(nsing1)- 1._r8)/alfmsh(nsing1) &
!              + widmsh(nsing1+1)*(alfmsh(nsing1+1)-1._r8)/alfmsh(nsing1+1)
!            print *,'zts(nsing1+1)',zts(nsing1+1)
             zts(nsing1+1) =  1._r8 
             zg(nsing1+1) = b*widmsh(nsing1+1)/alfmsh(nsing1+1)
             zd(nsing1+1) = zsing(nsing1+1)
             zc(nsing1+1) = zc(nsing1) - b*(alfmsh(nsing1+1)-1)* &
                            widmsh(nsing1+1)/alfmsh(nsing1+1)
!      
!            This is the check on the last zone to the edge.
             z = zts(nsing1)-zsing(nsing1)/b+zsing(nsing1+1)/b &
              + widmsh(nsing1)*(alfmsh(nsing1)- 1._r8)/alfmsh(nsing1) &
              + widmsh(nsing1+1)*(alfmsh(nsing1+1)-1._r8)/alfmsh(nsing1+1)
!
!            This is the check for no packing at the edge.
!            z = b + zc(nsing1)
!
             if( abs(z- 1._r8 ) >    1.E-10_r8 ) then
             write(itty, * )  ' ***error in nodes: zc(nosing)+b = ; ',z
             write(outmod, * )' ***error in nodes: zc(nosing)+b = ; ',z
             end if
           case(1)

             zts(nsing1+1) =  1._r8 
             zd(nsing1+1) = zsing(nsing1+1)

         end select
 
         jflag = 0
         do 88 i=1,nsing1
         z1 = zts(i+1) - zts(i)
         z2 = widmsh(i) + widmsh(i+1)
         if( z1 < z2) then
            write(itty,*)' *** i/M distance between rational surfaces = ', &
                 z1, ' < widmsh(',i+1,') + widmsh(',i,') = ' , z2
            write(itty,*)' *** check msin and/or decrease widmsh '
            write(outmod,*)' *** i/M distance between rational surfaces = ', &
                 z1, ' < widmsh(',i+1,') + widmsh(',i,') = ' , z2
            write(outmod,*)' *** check msin and/or decrease widmsh '
         widmsh(i) = widmsh(i)/(z2/z1 +  0.1_r8 )
         widmsh(i+1) = widmsh(i+1)/(z2/z1 +  0.1_r8 )
         jflag = jflag + 1
         end if
  88     continue
!
         jflag2 = jflag2 + 1 
         if( jflag /= 0  .AND.  jflag2 <= jfmax ) go to 100
!
!
!        Construct psi values at nodes between axis, surfaces and edge.
!        Each such region has three zones, 1, 2, and 3, extending in order
!        from the inner boundary to the outer boundary if each region.
!        Zones 1 and 3 are packed with a power via alfmsh, while Zone 2
!        is an equally spaced grid between them.
!
!         scan over regions...
!
         i3 = 0
         msing(1) = 1
         pkptsum=0
         do 85 i = 1,nsing1
           i0 = i3 + 1
           i1 = i0 + int( widmsh(i) * real(mminus) )
           i3 = int( zts(i + 1) *real(mminus) )  + 1
           msing(i+1) = i3 + 1 + pkptsum
           if( i == nsing1 ) msing(i+1) = mtot + 1
           i2 = i3 - int( widmsh(i+1) *real(mminus))
           if( i ==  nsing1 ) i2 = i2-1 !DPB the last zone near the edge
                                        ! needs a shifted index to have
                                        ! the derivative dzx/dz match at
                                        ! the zone boundary (one grid
                                        ! point further) 01/10
           ! mzoned is for output. 
           mzoned(1,i) = i1 - i0 + 1
           if( i /=  1 ) mzoned(1,i) = i1 - i0 + 1 + packpts(i)
           mzoned(2,i) = i2 - i1
           mzoned(3,i) = i3 - i2
           if( i /= nsing1 ) mzoned(3,i) = i3 - i2 + packpts(i+1)
           lzone(i) = 3
           zm = real( i3-i0+1 )
           if( (i /= 1)  .AND.  (i /= nsing1) ) zm = zm +  1.0_r8 
           js = i0 - 1
           if( i == 1 ) js = js + 1
           zdxs = zts(i+1) - zts(i)
!
!        construct zone 1
!
           do 86 in = i0,i1
             in2 = in + pkptsum
             z = zts(i) + zdxs*( real(in-js)/zm )
             zx = +zg(i  )*abs( (z-zts(i  ))/widmsh(i  ) )**alfmsh(i  ) &
                  + zd(i  )
             psinod(in2) = zx *psidif + psia(1)

86         continue
!
!        construct zone 2
!
           do 87 in = i1+1,i2
             in2 = in + pkptsum
             z = zts(i) + zdxs*( real(in-js)/zm )
             zx =  b*z + zc(i)
           psinod(in2) = zx *psidif + psia(1)
 
 87        continue
!
!        contruct zone 3
!  
           do 89 in = i2+1,i3
             in2 = in + pkptsum
             z = zts(i) + zdxs*( real(in-js)/zm )
!            only activate the packmth for the edge zone with packmth==1.
             if (packmth==1 .and. i/=nsing1) then
               pkmth=0
             else  
               pkmth=packmth
             endif
             select case(pkmth)
               case(0)
                 zx=-zg(i+1)*abs((z-zts(i+1))/widmsh(i+1))**alfmsh(i+1) &
                    + zd(i+1)
               case(1)
                 zx=c0-((1-af)*(-b*z/(c0+wx-1)**af-c1))**(1/(1-af));
             end select
 
            psinod(in2) = zx * psidif + psia(1)
 
 89       continue

          pkptsum=pkptsum+2*packpts(i) ! The number of surf packing points
 
 85       continue
 
!----------------------------------------------------------------
!
       else ! this is -M2 the hypergeometric distribution of grid points

         ! first set up the node density in each zone by chosing zts.
         ! Initially, set zts to zsing?
         zts(1) =  0._r8 
         do i=2,nsing1
            !zts(i) = zsing(i)**(3)
            !zts(i) = (real(i-1)/(nsing1))**(0.8)
            ! in this method the widmsh vector is used as input to specify
            ! the z locations of the surface zone delimitations.  This 
            ! allows for hand tuning the mesh to convergence.
            zts(i) = widmsh(i)
         end do
         zts(nsing1+1) =  1._r8 

         i3 = 0
         msing(1) = 1
         pkptsum=0
         do 91 i = 1,nsing1
           ! only i0 and i3 are needed, there is only one zone between surfs.
           i0 = i3 + 1
           i3 = int( zts(i + 1) *real(mminus) )  + 1
           msing(i+1) = i3 + 1 + pkptsum
           if( i == nsing1 ) msing(i+1) = mtot + 1

           ! mzoned is for output. report total in first zone
           mzoned(1,i) = i3 - i0 + 1
           mzoned(2,i) = 0
           mzoned(3,i) = 0
           !lzone(i) = 3
           zm = real( i3-i0+1 )
           if( (i /= 1) .AND. (i /= nsing1) ) zm = zm +  1.0_r8 
           js = i0 - 1
           !if( i == 1 ) js = js + 1
           zdxs = zts(i+1) - zts(i)

!          set up the normalization factor for the hypgeo function
           normfac=1.0_dp/beta(alfmsh(i),alfmsh(i+1))

           allocate(xhg(mzoned(1,i)))
           allocate(dz(mzoned(1,i)))
           allocate(z1hg(mzoned(1,i)))
           allocate(z2hg(mzoned(1,i)))
           allocate(z3hg(mzoned(1,i)))
           allocate(z4hg(mzoned(1,i)))
           allocate(z5hg(mzoned(1,i)))
           allocate(z6hg(mzoned(1,i)))
            
           if (i==1) then
             shft=1
           else
             shft=0
           end if
           do 92 in = i0,i3
             in2 = in + pkptsum
             xhg(in-js) = real(in-js-shft)/zm
             z = zts(i) + zdxs*( xhg(in-js) )

             dz(in-js)=(1./((zts(i+1)-zts(i))*mminus))* &
                       normfac*xhg(in-js)**(alfmsh(i)-1.0_dp)*  &
                       (1-xhg(in-js))**(alfmsh(i+1)-1.0_dp)
             z1hg(in-js)=hypgeo(cmplx(alfmsh(i),0.0_dp,kind=dpc), &
                                cmplx(1.0_dp-alfmsh(i+1),0.0_dp,kind=dpc), &
                                cmplx(1.0_dp+alfmsh(i),0.0_dp,kind=dpc), &
                              xhg(in-js))
             z3hg(in-js)=xhg(in-js)**(alfmsh(i))/alfmsh(i) 
             z4hg(in-js)=z1hg(in-js)*z3hg(in-js) 
             z5hg(in-js)=normfac*z4hg(in-js) 
             if (in-js>1) then
              z6hg(in-js)=(z5hg(in-js)-z5hg(in-js-1))*zm 
             endif
!            normfac*x(i)**(1+a)/(1+a)*hypgeo(1+a,-b,1+c,x(i)) 
 
             if (i==nsing1) then
               zx = z5hg(in-js)*(1.0_dp-zsing(i)) + zsing(i)
               dz(in-js)=dz(in-js)*(1.0_dp-zsing(i))
             else 
               zx = z5hg(in-js)*(zsing(i+1)-zsing(i)) + zsing(i)
               dz(in-js)=dz(in-js)*(zsing(i+1)-zsing(i))
             endif
             if (in2==1) then
               ! this has caused problems on some compilers (pgf90)
               ! why is not clear.  It's (0,0)^alfmsh =(0,0).  Just
               ! set it.
               ! write(*,*) zx,z3hg(in-js),xhg(in-js),alfmsh(i), &
               !           (0.0_r8,0.0_r8)**2.0_r8
               z3hg(in-js)=0.0_r8
               z4hg(in-js)=0.0_r8
               z5hg(in-js)=0.0_r8
               z6hg(in-js)=0.0_r8
               zx=0.0_r8 
             endif
             psinod(in2) = zx *psidif + psia(1)

!            write(98,*) z,real(dz(in-js)), &
!                        real(z1hg(in-js)),real(z3hg(in-js)), &
!                        real(z4hg(in-js)),real(z5hg(in-js)), &
!                        real(z6hg(in-js)),zx,psinod(in2)
 92        continue

           deallocate(xhg)
           deallocate(dz)
           deallocate(z1hg)
           deallocate(z2hg)
           deallocate(z3hg)
           deallocate(z4hg)
           deallocate(z5hg)
           deallocate(z6hg)

           pkptsum=pkptsum+2*packpts(i) ! The number of surf packing points
 91      continue
       end if
 101   continue

       if( nosing /= 0 ) then
!
! add equidistant surfaces on either side of rational surface.  The number
! of such points is controlled by packpts(i).  The points are laid out
! a distance dpsi=delpsa/2**(packpts(i)-j) from the rational surface on
! both sides, where j=1,packpts(i).  Thus, the total number within this
! group is 2*packpts(i).  For packpts(i)=1, the result is equivalent to
! the original method of having two points around the surface at a distance
! of dlay.  
!
      do 90 i = 1,nosing
      delpsl = psisin( i+1 ) - psinod( msing(i+1) - 1 )
      delpsr = psinod( msing(i+1) + 2*packpts(i) ) - psisin( i+1 )
      delpsi = delpsl/ 2._r8 
      if( delpsr < delpsl ) delpsi = delpsr/ 2._r8 
      delpsa = abs( dlay * psidif )
      if (dlaymth==1) then
        delpsi=delpsi/2._r8 ! in this method, we decrease again, but the
                            ! points remain adjusted relative to the local
                            ! grid spacing.
      else
        if( delpsa  <=  delpsi ) then
          delpsi = delpsa
        else
          if( imode == 1 ) then
            write(itty,  *)' ***nodes: delpsa  >   ' &
               ,'delpsi between psisin(',i,') and psisin(',i+1,')'
            write(itty,  *)' should decrease dlay'
            write(outmod,  *)' ***nodes: delpsa  >   ' &
               ,'delpsi between psisin(',i,') and psisin(',i+1,')'
            write(outmod,  *)' should decrease dlay'
          end if
        end if
      end if
!
!     Here we set the surface packing points as described above. 
!
!     first create a gaussian dzx using the psinod locations as a work space
!     and sumfac as a norm factor.
!
      sumfac=0
      do j=1,packpts(i)
        psinod(msing(i+1)-j  +packpts(i))= &
                              exp((real(j)/packpts(i))**2)
        sumfac=sumfac+psinod(msing(i+1)-j  +packpts(i))
      enddo
!     then divide the values by the (zx) norm factor 
      do j=1,packpts(i)
        psinod(msing(i+1)-j+packpts(i))= & 
                   psinod(msing(i+1)-j+packpts(i))/sumfac
      enddo
!     now crudely sum them up as an integration for the zx (being highly
!     accurate for zx is not crucial in this part of the algorithm) and
!     thus get the zx values in the space where delpsi=>1.  Thus here
!     the point at j=packpts(i) should be =1.
      if (packpts(i)>1) then
        do j=2,packpts(i)
          psinod(msing(i+1)-j+packpts(i))= &
             psinod(msing(i+1)-j    +packpts(i))+ &
             psinod(msing(i+1)-(j-1)+packpts(i))
        enddo
      endif
!     Now the first step dx is about the same size as the second.  But,
!     since there is no node placed AT the rational surface, we want 
!     the step size between the surface and the first point to be <=
!     half the next, so that the points surrounding the surface are
!     about equally spaced.  So, to do this we can simply shift the zx 
!     locations down by an operator zx=zx-(1-zx)*zx(1)/2 where zx(1)
!     is the first point, and the space is normalized as zx=>1 at delpsi 
!     so the point at delpsi is unchanged.
!
!     need to cycle down to get it right.
!
      do j=packpts(i),1,-1
         psinod(msing(i+1)-j+packpts(i))= &
           psinod(msing(i+1)-j+packpts(i))- &
             (1.0_r8-psinod(msing(i+1)-j+packpts(i)))* &
             psinod(msing(i+1)-1+packpts(i))/2
      enddo
!     Now the distribution is Gaussian-like-symmetric around the surface, 
!     with no point at the surface, and nearly equal spacing the vicinity 
!     of the surface.
!     load both sides of the surface from the work space
      do j=1,packpts(i)
        psinod(msing(i+1)+j-1+packpts(i))= psisin(i+1) + &
           delpsi*psinod(msing(i+1)-j+packpts(i))
        psinod(msing(i+1)-j  +packpts(i))= psisin(i+1) - &
           delpsi*psinod(msing(i+1)-j+packpts(i))
      enddo
       
     
!     sumfac=0
!     do j=1,packpts(i)
!       sumfac=sumfac+j
!     enddo
!     do j=1,packpts(i)
!       psinod(msing(i+1)-j  +packpts(i))= &
!                             psisin(i+1) - delpsi/2**(packpts(i)-j)
!       psinod(msing(i+1)+j-1+packpts(i))= &
!                             psisin(i+1) + delpsi/2**(packpts(i)-j)
!       jfac=0
!       do k=1,j
!         jfac=jfac+k
!       enddo
!       psinod(msing(i+1)-j  +packpts(i))= &
!                             psisin(i+1) - delpsi*jfac/sumfac
!       psinod(msing(i+1)+j-1+packpts(i))= &
!                             psisin(i+1) + delpsi*jfac/sumfac
!       psinod(msing(i+1)-j  +packpts(i))= &
!                             psisin(i+1) - j**2*delpsi/packpts(i)**2
!       psinod(msing(i+1)+j-1+packpts(i))= &
!                             psisin(i+1) + j**2*delpsi/packpts(i)**2
!     enddo

       if( imode == 1 ) then
      write(itty,8757)   delpsi 
      write(outmod,8757) delpsi 
 8757 format( ' half width of element near rat. surf = ',f13.8/)
       end if
!
 90   continue
!
      end if
!
 
!
!     construct new grid of psi surfaces for remapping onto later.
!
       msub = 2
!
!       reset m to reflect the new number of surfaces
!
       m = mtot
       mp = m + 1

!
       kn = 1
       km = 1
       psinew(1) = psinod(1)
        do 65 i = 1, mtot
       kn = kn + 1
       dps = ( psinod(kn)-psinod(kn-1) ) / msub(i)
       dtent(kn-1) = dps * msub(i)
       do 63 ik = 1, msub(i)
       km = km + 1
       psinew(km) = psinod(kn-1) + ik * dps
 63    continue
 65    continue
!
!      set nusurf
!
       nusurf = km
!
!      for diagnostics
!
      write(outmod,9007)(psisin(i),i=1,nsing1+1)
      write(outmod,9015)mtot
!
      write(itty  ,9017)
      write(outmod,9017)
      do i = 1, nsing1
      write(itty  ,9018) i, mzoned(1,i), mzoned(2,i), mzoned(3,i)
      write(outmod,9018) i, mzoned(1,i), mzoned(2,i), mzoned(3,i)
      end do	  
!
! set up table for plotting mesh characteristics
!
!      zl =  0._r8 
!      zr =  0._r8 
!      do 51 i = 1,mtot + 1
!      write(outnod,9020) i, psinod(i), zl, zr
!      zl =  2.0_r8  * abs( real((i+1)/2) - real(i+1)/ 2._r8  )
!      zr =  2.0_r8  * abs( real(i/2) - real(i)/ 2._r8  )
! 51   continue
!
      write(outmod,9010)(psinew(i),i=1,nusurf)
 9007 format(/,10(' psisin = ',7e12.5))
 9008 format(2x,'node',4x,'psinod',5x,'#1',2x,'#2')
 9020 format(2x,i4,e14.5,f4.0,f4.0)
 9010 format(' psinew = ',8e12.5)
!
!
 9012 format('this region is divided into ',i3,' zones',/ &
 ,' zone no.',   2x,' zone width',8x,'in psi',6x,' fin. elements',2x,' dpsi',/)
 9013 format(4x,i3,4x,f10.3,3x,e15.5,9x,i3,4x,e12.5)
 9014 format(/,' region no ',i3,' from ',e15.5,' to ',e15.5,/)
 9015 format(/,' there are', i5, ' finite elements in the mesh',/)
!
 9016 format('error in distributing the finite elements, check mzoned'    ,/ &
,' mtot = ',i5,' mzoned= ',5(11i4))
!
 9017 format(//' number of mesh nodes in:'/ &
         ' section     zone-1     zone-2     zone-3')
 9018 format(3x,i2,3x,3(8x,i3))
!
       if ( check1 ) then
          print *,'psia'
          print *,psia
          print *,'psinod'
          print *,psinod
          print *,'psinew'
          print *,psinew
       end if
          
!
9005     format(2x,' surf    psi-old      psi-new',/, &
     130(1x,i5,2e14.5,/) )
!
! check if psinew is in increasing order...
!
      jtest = 0
      do j = 1, nusurf-1
      if( psinew(j+1) <= psinew(j) ) then
      jtest = jtest+1
      write(  itty,*)' *** nodes: psinew not in increasing order' &
  , ' about node ',j,' psinew = ',psinew(j)
      write(outmod,*)' *** nodes: psinew not in increasing order' &
  , ' about node ',j,' psinew = ',psinew(j)
      write(  itty,*)'will stop '
      write(outmod,*)'will stop '
      write(*,*) 'nodes stop ',j,psinew(j+1),psinew(j)
      stop ' ERROR: psinew not in increasing order in nodes'
      end if
      end do
!
       return
       end subroutine


      subroutine resid_b(ivec,iok,x,ansr,ivecd,zinput,ninput,zoutput, &
                         noutput)
         implicit none
         integer ivec
         logical iok(ivec)  ! if iok(i) TRUE then SKIP element i of vector
         real*8 x(ivec)     ! input vector {x(i)}
         real*8 ansr(ivec)  ! output vector {f(x(i))}
         integer ivecd      ! auxilliary information vector size (.ge.ivec)
         integer ninput,noutput
         real*8 zinput(ivecd,ninput),zoutput(ivecd,noutput)

         real*8 cf,bs,we,cs,af,wx,c0

         bs=x(1)
         cf=zinput(1,1)
         we=zinput(1,2)
         cs=zinput(1,3)
         af=zinput(1,4)

         wx=1-bs*cs
         c0=(1+cf**(1/af)*(wx-1))/(1-cf**(1/af));

         ansr=-(c0-1+wx)**(1-af)/(1-af) + bs*we/(c0-1+wx)**af +  &
               (c0-1)**(1-af)/(1-af)

      end subroutine 

