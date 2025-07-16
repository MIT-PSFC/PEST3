subroutine eqi_geq_axis(nh,rgrid,nv,zgrid,psirz,nb,rbdy,zbdy, &
     rmaxis,zmaxis,Psiax,Psibdy,ierr)

  ! use Newton method to find magnetic axis at center min or max of Psi(R,Z)
  ! this is a root finder method using grad(Psi) -> 0 on axis

  ! return the adjusted axis location Psi(rmaxis,zmaxis) value as Psiax
  ! return the edge Psi value evaluate and averaged over the boundary contour

  !------------------
  ! mod DMC: Jan 2011: adding recovery code
  !   (deal with cases where 2d Newton search fails)
  !
  ! to do this, a shared module is created to contain Psi(R,Z) spline
  ! WARNING: thread safety not preserved, concurrent calls in a shared
  ! memory context not expected to survive...

  use eqi_geq_axis_data

  implicit NONE

  !------------------

  integer, intent(in) :: nh   ! R grid size & 1st dim of psiRZ
  integer, intent(in) :: nv   ! Z grid size & 2st dim of psiRZ

  real*8, intent(in) :: rgrid(nh)  ! R grid
  real*8, intent(in) :: zgrid(nv)  ! Z grid
  real*8, intent(in) :: psirz(nh,nv)  ! Psi(R,Z)

  integer, intent(in) :: nb   ! size of boundary contour
  real*8, intent(in) :: rbdy(nb),zbdy(nb)  ! the boundary contour

  !  on input, initial guess; on output, refined values:
  real*8, intent(inout) :: rmaxis  ! R of magnetic axis
  real*8, intent(inout) :: zmaxis  ! Z of magnetic axis

  real*8, intent(out) :: Psiax     ! Psi on axis (as adjusted)
  real*8, intent(out) :: Psibdy    ! Psi averaged along boundary

  integer, intent(out) :: ierr  ! status code 0=OK

  !------------------------
  ! local...
  !  spline information  -- arrays moved to module eqi_geq_axis_data

  integer :: ilinx,iliny,ict(6),idum

  !  spline lookup variables (for bdy eval)

  integer, dimension(:), allocatable :: iia,jja
  real*8, dimension(:), allocatable :: rparama,zparama,hra,hria,hza,hzia
  real*8, dimension(:), allocatable :: psiaa,raa,zaa,wk

  ! Newton search variables
  integer :: iter,itmax,i,j,imin,jmin
  real*8 :: wkrmag(1),wkzmag(1),psidat(1,6)
  real*8 :: ztol,dpsidR,dpsidZ,d2psidR2,d2psidZ2,zrmp,zzmp,zdist,zmin_psi
  real*8 :: zdpsi_sum
  logical :: iexit

  integer :: nbdy,ib,jb

  real*8 xdum,x1(1),x2(1),zinput(1,2),zoutput(1,1)

  real*8 :: Rmin,Rmax,Zmin,Zmax

  logical :: iok0(1)

  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(12,100)
  real*8, parameter :: ZERO = 0.0_R8

  EXTERNAL eqi_geq_axis_Rsrch

  !-------------------------
  ! Newton iteration to adjust magnetic axis to match the bicubic
  ! spline representation over the user specified grids
  !   (flag indicate these differ from the xplasma Psi(R,Z) internal
  !   grids)

#ifdef __DEBUG
  allocate(wk(nv)); wk=0.0d0

  call r8_grf3f1(rgrid,psirz,zgrid,wk,nh,nv,nh, &
       'R','Psi(R,Z)','Z','m','Wb/rad','m', &
       'eqi_geq_axis.f90 debug plot','eqi_geq_axis Psi(R,Z)',0)

  deallocate(wk)
#endif

  Psiax = ZERO
  Psibdy = ZERO
  ierr = 0

  xdum = ZERO

  Rmin = Rgrid(1)
  Rmax = Rgrid(nh)

  Zmin = Zgrid(1)
  Zmax = Zgrid(nv)

  zdist = min((Rmax-Rmin),(Zmax-Zmin))/10

  ztol = 1.0d-12*max(Rmax,max(abs(Zmin),abs(Zmax)))

  ! set up bicubic spline and (R,Z) lookup data

  allocate(rpkg(nh,4),zpkg(nv,4),psirz_spl(4,nh,nv))

  call r8genxpkg(nh,rgrid,rpkg,0,0,0,xdum,3,ierr)
  if(ierr.ne.0) go to 990

  call r8genxpkg(nv,zgrid,zpkg,0,0,0,xdum,3,ierr)
  if(ierr.ne.0) go to 990

  psirz_spl(1,:,:) = psirz(:,:)
  call r8mkbicub(rgrid,nh,zgrid,nv,psirz_spl,nh, &
       0,xdum,0,xdum, 0,xdum,0,xdum, &
       ilinx,iliny, ierr)
  if(ierr.ne.0) go to 990

  ! set up for evaluation of Newton iteration data
  ict = (/ 1,1,1,1,1,0 /)  ! value & 1st & 2nd derivatives
  ! (cross derivative omitted as unnecesary in this context)

  psidat = ZERO
  wkrmag = rmaxis
  wkzmag = zmaxis

  idum = 0

  iexit = .FALSE.
  iter = 0
  itmax = 100  ! max iteration count

  do
     iter = iter + 1
     if(iter.gt.itmax) then
        ierr = 1
#ifdef __DEBUG
        write(6,*) ' *** eqi_geq_axis: 2d Newton loop did not converge *** '
#endif
        exit
     endif

     ! save old value
     zrmp=wkrmag(1)
     zzmp=wkzmag(1)

     call r8xlookup(1,wkrmag,nh,rpkg,2,ii,rparam,hr,hri,idum)
     call r8xlookup(1,wkzmag,nv,zpkg,2,jj,zparam,hz,hzi,idum)

     ! evaluate bicubic spline representation of Psi and its derivatives

     call r8fvbicub(ict,1,1, &
          psidat,ii,jj,rparam,zparam,hr,hri,hz,hzi, &
          psirz_spl,nh,nv)

     ! extract value
     Psiax = psidat(1,1)

     ! extract derivatives
     dPsidR = psidat(1,2)
     dPsidZ = psidat(1,3)
     d2PsidR2 = psidat(1,4)
     d2PsidZ2 = psidat(1,5)

     ! if converged, exit now
     if(iexit) exit

     ! done if 1st derivatives are both zero
     if((dPsidR.eq.ZERO).AND.(dPsidZ.eq.ZERO)) exit

     ! Newton step in R direction
     !   error if 2nd derivative is zero

     if(d2PsidR2.eq.ZERO) then
        iter=itmax
        cycle
     endif

     ! new Rmag
     if(dPsidR.ne.ZERO) then
        wkrmag(1) = zrmp - dPsidR/d2PsidR2
     endif

     ! Newton step in Z direction
     !   error if 2nd derivative is zero

     if(d2PsidZ2.eq.ZERO) then
        iter=itmax
        cycle
     endif

     ! new Zmag
     if(dPsidZ.ne.ZERO) then
        wkzmag(1) = zzmp - dPsidZ/d2PsidZ2
     endif

     ! compare changes to tolerance
     iexit = (max(abs(wkzmag(1)-zzmp),abs(wkrmag(1)-zrmp)).lt.ztol)

     ! force early exit if wandered outside the grid

     if(wkRmag(1).lt.Rmin) iter=itmax
     if(wkRmag(1).gt.Rmax) iter=itmax
     if(wkZmag(1).lt.Zmin) iter=itmax
     if(wkZmag(1).gt.Zmax) iter=itmax

  enddo

  if(ierr.ne.0) then
     ! Newton loop not converged...
     ! try Ridder search

     eqi_eps = 0.000001d0*min((Rgrid(2)-Rgrid(1)),(Zgrid(2)-Zgrid(1)))
     zmin_psi = 1.0d30
     do j=1,nv
        do i=1,nh
           if((abs(Rgrid(i)-Rmaxis).lt.zdist).AND. &
                (abs(Zgrid(j)-Zmaxis).lt.zdist)) then
              if(psiRZ(i,j).lt.zmin_psi) then
                 zmin_psi = psiRZ(i,j)
                 imin=i
                 jmin=j
              endif
           endif
        enddo
     enddo

     zdpsi_sum = &
          (psiRZ(imin+1,jmin)-psiRZ(imin,jmin))/(Rgrid(imin+1)-Rgrid(imin)) +&
          (psiRZ(imin,jmin)-psiRZ(imin-1,jmin))/(Rgrid(imin)-Rgrid(imin-1)) +&
          (psiRZ(imin,jmin+1)-psiRZ(imin,jmin))/(Zgrid(jmin+1)-Zgrid(jmin)) +&
          (psiRZ(imin,jmin)-psiRZ(imin,jmin-1))/(Zgrid(jmin)-Zgrid(jmin-1))

     zdpsi_sum = zdpsi_sum/4
     eqi_eta = 0.001d0*zdpsi_sum

     iok0(1) = .FALSE.
     x1(1) = Rgrid(imin-3)
     x2(1) = Rgrid(imin+3)
     zinput(1,1)=Zgrid(jmin-3)
     zinput(1,2)=Zgrid(jmin+3)

     call zriddery(1,iok0,x1,x2,eqi_eps,eqi_eta, &
          eqi_geq_axis_Rsrch,wkRmag,ierr, &
          1,zinput,2,zoutput,1)

     if(ierr.ne.0) then
#ifdef __DEBUG
        write(6,*) ' *** eqi_geq_axis: 2d Ridder search failed *** '
#endif
        go to 990
     else
#ifdef __DEBUG
        write(6,*) ' but eqi_geq_axis: 2d Ridder search OK. '
#endif
        continue
     endif

     wkZmag = zoutput(1,1)

     call r8xlookup(1,wkrmag,nh,rpkg,2,ii,rparam,hr,hri,idum)
     call r8xlookup(1,wkzmag,nv,zpkg,2,jj,zparam,hz,hzi,idum)

     ! evaluate bicubic spline representation of Psi and its derivatives

     call r8fvbicub(ict,1,1, &
          psidat,ii,jj,rparam,zparam,hr,hri,hz,hzi, &
          psirz_spl,nh,nv)

     ! extract value
     Psiax = psidat(1,1)

  endif

  ! normal exit from Newton loop;
  ! converged: new mag axis for this representation

  rmaxis = wkRmag(1)
  zmaxis = wkZmag(1)

  ! Psiax was evaluated in the Newton loop.

  ! Now evaluate Psibdy by averaging over the boundary contour

  nbdy = nb
  if((rbdy(1).eq.rbdy(nb)).AND.(zbdy(1).eq.zbdy(nb))) nbdy = nbdy - 1

  nbdy = nbdy * 2

  allocate(iia(nbdy),jja(nbdy),rparama(nbdy),zparama(nbdy))
  allocate(hra(nbdy),hria(nbdy),hza(nbdy),hzia(nbdy))

  allocate(psiaa(nbdy),raa(nbdy),zaa(nbdy))
  
  jb=0
  do ib=1,nbdy/2
     jb=jb + 1
     raa(jb) = rbdy(ib)
     zaa(jb) = zbdy(ib)
     jb=jb + 1
     if(ib.lt.nb) then
        raa(jb) = (rbdy(ib)+rbdy(ib+1))/2
        zaa(jb) = (zbdy(ib)+zbdy(ib+1))/2
     else
        raa(jb) = (rbdy(ib)+rbdy(1))/2
        zaa(jb) = (zbdy(ib)+zbdy(1))/2
     endif
  enddo

  call r8xlookup(nbdy,raa,nh,rpkg,2,iia,rparama,hra,hria,idum)
  call r8xlookup(nbdy,zaa,nv,zpkg,2,jja,zparama,hza,hzia,idum)

  ! set up for evaluation of Psi values only
  ict = (/ 1,0,0,0,0,0 /)

  call r8fvbicub(ict,nbdy,nbdy, &
       psiaa,iia,jja,rparama,zparama,hra,hria,hza,hzia, &
       psirz_spl,nh,nv)

  Psibdy = ZERO
  do ib=1,nbdy
     Psibdy = Psibdy + psiaa(ib)
  enddo
  Psibdy = Psibdy/nbdy

  deallocate(iia,jja,rparama,zparama,hra,hria,hza,hzia)
  deallocate(psiaa,raa,zaa)

990 continue
  deallocate(rpkg,zpkg,psirz_spl)

end subroutine eqi_geq_axis
