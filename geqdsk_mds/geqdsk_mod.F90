module geqdsk_mod
!
! package geqdsk data
! pletzer@ pppl.gov Mon Jul 31 16:04:04 EDT 2000
! dmccune@pppl.gov -- Jan 2001 -- added MDSplus EFIT tree access
!
! Refer to http://lithos.gat.com/efit/g_eqdsk.html for a detailed
! description of the data contained in g-eqdsk files. The module
! herein follows the terminology described in the above document; variables
! bear the same name with their units trailing after a underscore (_). A
! double underscore (__) stands for /, as in Wb__Rad (=Weber/rads).
 
  implicit NONE

  integer, parameter :: geq_r8 = selected_real_kind(12,100)
  real(geq_r8), parameter :: geq_pi = 3.141592653589793116
  real(geq_r8), parameter :: geq_2pi = 6.2831853071795862320
  real(geq_r8), parameter :: geq_mu0 = 1.2566370614359172885e-06
 
  type geqdsk
 
 
     integer  :: NW ! Number of horizontal R grid points
     integer  :: NH ! Number of vertical Z grid points
     character*8 :: CASE_(6)
     real(geq_r8) :: RDIM_M, ZDIM_M ! box dimensions
     real(geq_r8) :: RCENTR_M ! Major radius (not magnetic axis!)
     real(geq_r8) :: RLEFT_M ! left side of computational box
     real(geq_r8) :: ZMID_M ! mid elevation of computational box
     real(geq_r8) :: RMAXIS_M, ZMAXIS_M ! Magnetic axis position
     real(geq_r8) :: BCENTR_T, CURRENT_A
     ! poloidal flux/2*pi on magnetic axis and plasma boundary
     real(geq_r8) :: SIMAG_Wb__Rad, SIBRY_Wb__Rad
     real(geq_r8), dimension(:), allocatable :: FPOL_TM ! toroidal current
     real(geq_r8), dimension(:), allocatable :: PRES_Nt__M2 ! pressure
     real(geq_r8), dimension(:), allocatable :: FFPRIM_T2M2Rad__Wb ! F F'(psi)
     real(geq_r8), dimension(:), allocatable :: pprime_NtRad__M2Wb ! p'(psi)
     real(geq_r8), dimension(:,:), allocatable :: psirz_Wb__Rad
     real(geq_r8), dimension(:), allocatable :: QPSI
     integer  :: nbbbs ! Number of plasma-vacuum boundary points
     integer  :: limitr ! Number of limiter points
     real(geq_r8), dimension(:), allocatable :: RBBBS_M, ZBBBS_M ! Plasma-vacuum
     real(geq_r8), dimension(:), allocatable :: RLIM_M, ZLIM_M ! Limiter geometry
 
!  to support MDSplus... (MDSOPT in reserve, not currently used)
 
     integer :: access_method ! 1:file, 2:MDS+ unreduced, 3: MDS+ reduced
     character*180 :: IDENT ! data source (filename or MDSplus path)
     real(geq_r8) :: MDSTIME ! current time of interest (if MDSplus)
     character*10 :: MDSOPT ! auxilliary MDSplus data processing option
     character*10 :: MDSTOK ! experiment (tokamak) machine name
     integer MDS_SHOT        ! shot number
     logical MDSLICE         ! .TRUE. to use slice-by-slice MDSplus access
 
     character*10 :: marker  ! set to "xplasma:ok" if setup is successful.
 
!
!  the following are used when accessing an MDSplus TREE.  These
!  arrays cache the time dependent EFIT data; only "valid" timepoints
!  are retained.  see the contained routines.  DMC Jan 2001.
!
 
     integer :: T_NUM_VALID
     real(geq_r8), dimension(:), allocatable :: GTIME_S
                     ! (if MDSplus EFIT tree):  "actual" timebase, seconds
 
     integer, dimension(:), allocatable :: T_INDEX  ! into "original" timebase

     real(geq_r8) :: CURTIME  ! time (s) of MDS+ slice currently selected
 
!  the "original" timebase contains some points for which there is no valid
!    data, and these are excluded from the "actual" timebase.  T_INDEX
!    indexes from the "actual" to the "orginal" time dimension.
 
     integer, dimension(:), allocatable :: T_NW
     integer, dimension(:), allocatable :: T_NH
 
     integer :: kt_max  ! max index for T_RDIM,T_ZDIM ....
 
     integer :: kt_sw2  ! storage ordering for profiles f(psi,t)
     integer :: kt_sw3  ! storage ordering for profile psi(R,Z,t)
 
! kt_max:
! (for some EFIT trees, dimensioning scalars are time dependent
!     => kt_max=t_num_valid;
!  for others they are not => kt_max=1)
!
! affected quantities:  T_NW,T_NH,T_RDIM,T_ZDIM,T_RCENTR,T_RLEFT,T_CMID
!
! kt_sw2,kt_sw3:
! (for some EFIT trees, profiles are stored psi-contiguous; and in
! others, t-contiguous; the code has to adapt)
 
     real(geq_r8), dimension(:), allocatable :: T_RDIM,T_ZDIM
     real(geq_r8), dimension(:), allocatable :: T_RCENTR
     real(geq_r8), dimension(:), allocatable :: T_RLEFT
     real(geq_r8), dimension(:), allocatable :: T_ZMID
     real(geq_r8), dimension(:), allocatable :: T_RMAXIS, T_ZMAXIS
     real(geq_r8), dimension(:), allocatable :: T_BCENTR, T_CURRENT
     real(geq_r8), dimension(:), allocatable :: T_SIMAG, T_SIBRY
 
     real(geq_r8), dimension(:,:), allocatable :: T_FPOL
     real(geq_r8), dimension(:,:), allocatable :: T_PRES
     real(geq_r8), dimension(:,:), allocatable :: T_FFPRIM
     real(geq_r8), dimension(:,:), allocatable :: T_PPRIME
     real(geq_r8), dimension(:,:), allocatable :: T_QPSI
 
     real(geq_r8), dimension(:,:,:), allocatable :: T_PSIRZ
 
     integer, dimension(:), allocatable :: T_NBBBS
     real(geq_r8), dimension(:,:), allocatable :: T_RBBBS,T_ZBBBS
 
!    the "time dependent" limiter data is not read; the limiter is assumed
!    not to move.
 
  end type geqdsk
 
! static items:
 
contains
 
  subroutine geq_xz_is_in_plasma(geq, x, z, iok, reltol)
 
    ! Return iok=0 if (x,z) is within plasma boundary
 
    ! Method: integrate over the plasma boundary,
    ! the angle between the point (x,z)
    ! and a segment of the boundary. If the result is
    ! nearly 2*pi then the point must be inside, if
    ! it's close to zero it must be outside.
 
    ! (dmc) also compute distance from (x,z) to the nearest point
    ! of approach to the boundary (treated as a sequence of line
    ! segments), normalized to the distance of the boundary point
    ! furthest from (x,z)
 
    ! iok = 0 => point inside
    ! iok = 1 => point ouside
    ! iok = 2 => can't tell for sure
 
    ! reltol
 
    implicit none
    type(geqdsk), intent(in) :: geq
    real(geq_r8), intent(in) :: x, z
    integer, intent(out) :: iok
    real(geq_r8), intent(out) :: reltol
 
    real(geq_r8) angle, dotprod, crsprod, x1, x2, z1, z2, tol
    real(geq_r8) distmin2,distmax2,dist2,dist12sq,alpha,beta
    real(geq_r8), parameter :: ZERO = 0.0d0
    real(geq_r8), parameter :: ONE  = 1.0d0
 
    integer i
 
    !------------------------------------------
 
    iok = 1
    reltol = 1000.0_geq_r8
 
    angle = 0._geq_r8
 
    x1=geq%RBBBS_M(geq%NBBBS)
    z1=geq%ZBBBS_M(geq%NBBBS)
 
    ! reference distance:  plasma bdy:  min((Rmax-Rmin),(Zmax-Zmin))
    distmax2=min( &
         maxval(geq%RBBBS_M(1:geq%NBBBS))-minval(geq%RBBBS_M(1:geq%NBBBS)), &
         maxval(geq%ZBBBS_M(1:geq%NBBBS))-minval(geq%ZBBBS_M(1:geq%NBBBS)) )
    distmax2=distmax2*distmax2  ! ref. dist. squared
 
    if(distmax2.eq.0._geq_r8) return   ! assume outside -- degenerate case.
 
    dist2=(x-x1)*(x-x1)+(z-z1)*(z-z1)
    distmin2=dist2
 
    do i = 1, geq%NBBBS - 1
 
       x1 = geq%RBBBS_M(i  ); z1 = geq%ZBBBS_M(i  )
       x2 = geq%RBBBS_M(i+1); z2 = geq%ZBBBS_M(i+1)
 
       dist2=(x-x1)*(x-x1)+(z-z1)*(z-z1)
       distmin2=min(distmin2,dist2)  ! find distance^2 to closest bdy pt
 
       dotprod = (x2-x)*(x1-x) + (z2-z)*(z1-z)
       crsProd = (x2-x)*(z1-z) - (z2-z)*(x1-x)
       dist12sq = (x2-x1)*(x2-x1)+(z2-z1)*(z2-z1)
       !  skip coincident points
       if((crsProd.eq.ZERO).or.(dist12sq.le.ZERO)) cycle
 
       angle = angle + abs(atan2(crsProd, dotprod))
 
       !  refine distance estimate by considering nearest approach to
       !  line segment itself, not just the end pts.
       !     P=(x,z)  P1=(x1,z1)  P2=(x2,z2)  DEL=(P2-P1)
       !     DELperp = vector perpendicular to DEL and of same length
       !
       !  bdy segment is represented by Pseg = P1 + alpha*DEL
       !     with alpha in [0,1]
       !
       !  line from P normal to Pseg is Plin = P + beta*DELperp
       !     and the intersection is found by setting Plin=Pseg and
       !     solving the 2x2 system for alpha & beta...
       !     if alpha in [0,1] then P is closer to the segment than
       !     to the segment end points, and the distance^2 is
       !     abs(beta)*|DEL|.
 
       alpha=-( (x2-x1)*(x1-x)+(z2-z1)*(z1-z) )/dist12sq
       beta = ( (x2-x1)*(z1-z)-(z2-z1)*(x1-x) )/dist12sq
 
       !debug  write(6,*) x1+alpha*(x2-x1), z1+alpha*(z2-z1)
       !debug  write(6,*) x -beta*(z2-z1),  z +beta*(x2-x1)
 
       if((alpha.ge.ZERO).and.(alpha.le.ONE)) then
          distmin2=min(distmin2,(dist12sq*beta*beta))
       endif
 
    enddo
 
 
    ! may need to adjust this tol
    tol = 0.01_geq_r8 * (2._geq_r8 * geq_pi)
 
    iok = 2
    if(abs(angle - geq_2pi) < tol) iok = 0
    if(abs(angle          ) < tol) iok = 1
 
    reltol=sqrt(distmin2/distmax2)
 
  end subroutine geq_xz_is_in_plasma
 
!------------------------------------------------------
! compute GS error
 
  subroutine geq_GSerror(isign, geq, nx1, nz1, x, z, bdist, gserror, ier)
 
    ! Compute rel. Grad Shafranov error *inside* plasma boundary
    ! (error is set to zero outside) (also, points "too close" to
    ! the boundary have error set to zero, see "bdist" control).
    !   ier=1 is set if no points are found (which generally indicates a
    !   very bad input dataset).
 
    use ezspline_obj
    use ezspline
    use geqdsk_aux

    implicit none
 
    integer, intent(in)      :: isign    ! +1 for EFIT G-eqdsk

    ! if Psi(a) < Psi(0), isign = -1 -- never true for EFIT but it
    !  is true for some imitators...

    type(geqdsk), intent(in) :: geq      ! G-eqdsk object
    integer, intent(in)      :: nx1, nz1 ! X,Z grids...
    real(geq_r8), intent(in) :: x(nx1), z(nz1)
 
    real(geq_r8), intent(in) :: bdist    ! min. relative distance to bdy
 
    real(geq_r8), intent(out) :: gserror(nx1, nz1)
    integer, intent(out) :: ier
 
    !  points within bdist*min((Rmax-Rmin),(Zmax-Zmin)) of the boundary
    !  have their error set to zero.  {Rmin,Rmax,Zmin,Zmax} are taken
    !  from the G-Eqdsk plasma boundary data.
 
    type(ezspline2_r8) :: psispl
    type(ezspline1_r8) :: ppspl, ggpspl
    real(geq_r8), dimension(:,:), allocatable :: d2xPsi, dxPsi, d2zPsi
    real(geq_r8), dimension(:), allocatable :: zwk
    real(geq_r8) :: pprime, ggprime, psi, r_jphi0
    real(geq_r8) :: dist_norm
    real(geq_r8) :: loc_psimin,loc_psimax
    integer iok, i, j, idum1, idum2, ict, ictbb
    character*30 zlabel       ! (for debugging)
 
    iok = 0
    ier = 0
 
    ! compute del^* psi
 
    call ezspline_init(psispl, geq%NW, geq%NH, (/0,0/), (/0,0/), iok)
    if(iok/=0) call ezspline_error(iok)
 
    psispl%x1 = (/ &
         & ( geq%RLEFT_M + &
         & geq%RDIM_M*real(i-1,geq_r8)/real(geq%NW-1, geq_r8), &
         & i = 1, geq%NW)/)
    psispl%x2 = (/ &
         & ( geq%ZMID_M - geq%ZDIM_M/2._geq_r8 + &
         & geq%ZDIM_M*real(i-1,geq_r8)/real(geq%NH-1, geq_r8), &
         & i = 1, geq%NH)/)
 
    call ezspline_setup(psispl, isign*geq%psirz_Wb__Rad, iok)
    if(iok/=0) call ezspline_error(iok)
 
    allocate(d2xPsi(nx1, nz1), stat=iok)
    allocate(dxPsi(nx1, nz1), stat=iok)
    allocate(d2zPsi(nx1, nz1), stat=iok)
 
    call ezspline_derivative(psispl, 1, 0, nx1, nz1, x, z, dxPsi, iok)
    if(iok/=0) call ezspline_error(iok)
 
    call ezspline_derivative(psispl, 2, 0, nx1, nz1, x, z, d2xPsi, iok)
    if(iok/=0) call ezspline_error(iok)
 
    call ezspline_derivative(psispl, 0, 2, nx1, nz1, x, z, d2zPsi, iok)
    if(iok/=0) call ezspline_error(iok)
 
    gserror = 0._geq_r8 ! d2xPsi - dxPsi/x + d2zPsi ! so far only del^* psi
 
    ! compute R^2 p' + g g'
 
    call ezspline_init(ppspl, geq%NW, (/0,0/), iok)
    if(iok/=0) call ezspline_error(iok)
 
    ppspl%x1 = (/ (  geq%SIMAG_Wb__Rad + (geq%SIBRY_Wb__Rad-geq%SIMAG_Wb__Rad) * &
         & real(i-1,geq_r8)/real(geq%NW-1), i = 1, geq%NW) /)
    ppspl%x1 = isign*ppspl%x1
 
    call ezspline_setup(ppspl, isign*geq_mu0*geq%pprime_NtRad__M2Wb, iok)
    if(iok/=0) call ezspline_error(iok)
 
    call ezspline_init(ggpspl, geq%NW, (/0,0/), iok)
    if(iok/=0) call ezspline_error(iok)
 
    ggpspl%x1 = ppspl%x1
 
    call ezspline_setup(ggpspl, isign*geq%FFPRIM_T2M2Rad__Wb, iok)
    if(iok/=0) call ezspline_error(iok)
 
    ! normalization = "nominal r_jphi"
    r_jphi0 = maxval(abs( &
         geq%RMAXIS_M**2 * geq_mu0*geq%pprime_NtRad__M2Wb(1:geq%NW) + &
         geq%FFPRIM_T2M2Rad__Wb(1:geq%NW) ))
    r_jphi0 = isign * r_jphi0
 
    ict = 0
    ictbb = 0
    do j = 1, nz1
       do i = 1, nx1
          ! get psi
          call ezspline_interp(psispl, x(i), z(j), psi, iok)

          if((i*j).eq.1) then
             loc_psimin=psi
             loc_psimax=psi
          else
             loc_psimin=min(loc_psimin,psi)
             loc_psimax=max(loc_psimax,psi)
          endif

          if( (psi>isign*geq%SIMAG_Wb__Rad .and. psi<isign*geq%SIBRY_Wb__Rad)  .or.  &
              (psi>isign*geq%SIBRY_Wb__Rad .and. psi<isign*geq%SIMAG_Wb__Rad) ) then
             ictbb = ictbb + 1

             ! most likely inside the plasma
             iok = 1 ! outside by default
             call geq_xz_is_in_plasma(geq, x(i), z(j), iok, dist_norm)
             if(iok==0 .and. dist_norm > bdist ) then
                ! definitely inside the plasma
                ict = ict + 1
                call ezspline_interp(ppspl, psi, pprime, iok)
                call ezspline_interp(ggpspl, psi, ggprime, iok)
                gserror(i,j) = (d2xPsi(i,j) - dxPsi(i,j)/x(i) + d2zPsi(i,j) &
                     & + x(i)**2 * pprime + ggprime)/r_jphi0
             endif
          endif
       enddo
    enddo
 
#ifdef __DEBUG

 
  !! Allocate(Zwk(max(nx1,nz1))); zwk=0
  !! zlabel='geqdsk_mds library (geqdsk_mod.f90)'
 
  !! call r8_grf3f1(x,gserror,z,zwk,nx1,nz1,nx1, &
  !!     'R','Normalized GS error','Z','m',' ','m', &
  !!     'eqi_fromgeqdsk debug plot',zlabel,0)
 
  !! deallocate(zwk)


#endif
    call ezspline_free(psispl, iok)
    if(iok/=0) call ezspline_error(iok)
    call ezspline_free(ppspl, iok)
    if(iok/=0) call ezspline_error(iok)
    call ezspline_free(ggpspl, iok)
    if(iok/=0) call ezspline_error(iok)
 
    deallocate(d2xPsi)
    deallocate(dxPsi)
    deallocate(d2zPsi)
 
    if(ict.eq.0) then
       ier=1
       if(isign.ne.1) then
          write(lunmsg,*) ' %geq_GSerror: EQDSK format variant detected.'
       endif

       if(ictbb.eq.0) then
          write(lunmsg,*) ' ?geq_GSerror: no Psi(R,Z) values in range:'
          write(lunmsg,*) '   (Psi values expected to be in Wb/rad)'
          write(lunmsg,*) '   Psi(axis) data: ',isign*geq%SIMAG_Wb__Rad
          write(lunmsg,*) '   Psi(bdy) data:  ',isign*geq%SIBRY_Wb__Rad
          write(lunmsg,*) '   min value in Psi(R,Z): ',loc_psimin
          write(lunmsg,*) '   max value in Psi(R,Z): ',loc_psimax
       else
          write(lunmsg,*) &
               ' ?geq_GSerror: valid Psi(R,Z) values not found inside bdy'
       endif
    endif

  end subroutine geq_GSerror
 
 !------------------------------------------------------
! set LUN for G-EQdisk file i/o
 
  subroutine geq_setlun(ilun)
 
    use geqdsk_aux
 
    implicit none
    integer, intent(in) :: ilun    ! unit number
 
    neqdsk=ilun
 
    return
  end subroutine geq_setlun
 
!------------------------------------------------------
! retrieve LUN for G-EQdisk file i/o
 
  subroutine geq_getlun(ilun)
 
    use geqdsk_aux
 
    implicit none
    integer, intent(out) :: ilun    ! unit number
 
    ilun=neqdsk
 
    return
  end subroutine geq_getlun
 
!------------------------------------------------------
  subroutine geq_writ(geq, filename, ier)
!
! ...write a G-EQDSK file from a stored objec...
!
! this could be used e.g. to extract files from remote MDSplus servers
!
    use geqdsk_aux
 
    implicit none
    type(geqdsk), intent(in) :: geq
    character*(*), intent(in) :: filename
    integer, intent(out) :: ier
 
! local stuff
 
    integer iok
    integer :: idum = 0
    real(geq_r8) :: adum = 0.0_geq_r8
    integer i,j
 
! executable code
 
    ier = 0
 
    if((geq%marker.ne.'xplasma:ok').AND.(geq%marker.ne.'scrunch2')) then
       write(lunmsg,*) ' %geq_writ: geqdsk object "geq" not ready.'
       write(lunmsg,*) '  attempt to write G-EQDSK file aborted.'
       ier = 1
       return
    endif
 
    open(neqdsk, file=filename, status='unknown', action='write', iostat=iok)
 
    if (iok /= 0) then
       write(lunmsg,*) ' %geq_writ: open failure: ',filename
       ier = 1
       return
    endif
 
    write (neqdsk,2000) (geq%case_(i),i=1,6),idum,geq%nw,geq%nh
    write (neqdsk,2020) geq%rdim_m, geq%zdim_m, &
         & geq%rcentr_m, geq%rleft_m, geq%zmid_m
    write (neqdsk,2020) geq%rmaxis_m,geq%zmaxis_m, &
         & geq%simag_Wb__Rad,geq%sibry_Wb__Rad,geq%bcentr_t
    write (neqdsk,2020) geq%current_a,geq%simag_Wb__Rad, &
         & adum,geq%rmaxis_m,adum
    write (neqdsk,2020) geq%zmaxis_m,adum,geq%sibry_Wb__Rad,adum,adum
 
    write (neqdsk,2020) (geq%fpol_tm(i),i=1,geq%nw)
    write (neqdsk,2020) (geq%pres_Nt__M2(i),i=1,geq%nw)
    write (neqdsk,2020) (geq%ffprim_T2M2Rad__Wb(i),i=1,geq%nw)
    write (neqdsk,2020) (geq%pprime_NtRad__M2Wb(i),i=1,geq%nw)
    write (neqdsk,2020) ((geq%psirz_Wb__Rad(i,j),i=1,geq%nw),j=1,geq%nh)
    write (neqdsk,2020) (geq%qpsi(i),i=1,geq%nw)
    write (neqdsk,2022) geq%nbbbs,geq%limitr
 
    write (neqdsk,2020) (geq%rbbbs_m(i),geq%zbbbs_m(i),i=1,geq%nbbbs)
 
    write (neqdsk,2020) (geq%rlim_m(i),geq%zlim_m(i),i=1,geq%limitr)
 
    close(neqdsk)
 
2000 format (6a8,3i4)
2020 format (5e16.9)
2022 format (2i5)
 
  end subroutine geq_writ
!------------------------------------------------------
! initialize a G-EQdisk object from file or MDS+ data
 
  subroutine geq_init(geq, filenami, ier)
 
    use geqdsk_aux
 
    implicit none
    type(geqdsk), intent(inout) :: geq
    character*(*), intent(in) :: filenami
    integer, intent(out) :: ier
 
! local stuff
 
    character*200 filename,zpath,ztail
    integer iok, i, j, idum, istat, indx, ilast
    character*4 mtest
    real(geq_r8) adum,zsimag
 
! local, for MDSplus support
 
    logical ireduce
    integer icln,ict,ilen,str_length,ishot,iopen
    character*10 opt
    character*11 itest
    character*50 mds_server
    character*20 mds_tree
    character*10 mds_id1
    character*10 mds_id2
    character*20 mds_tstr,z20
    real(geq_r8) ztime
 
    character*180 ident
    character*48 tmplbl
    integer iltmp,iltmp2
 
!-------------------------------------
! executable code
 
    ier = 0

! dmc modification -- translate EFIT_INDEX file specification if given
    itest=filenami(1:11)
    call uupper(itest)
    if(itest.eq.'EFIT_INDEX:') then
       call geqdsk_findex(filenami,ilast,indx)
       call geqdsk_index1_path(neqdsk,filenami(1:ilast),zpath)
       if(zpath.eq.'?') then
          write(lunmsg,*) ' ?geqdsk_mod.geq_init: index file access failure.'
          ier=1
       else
          call geqdsk_index1_gfile(neqdsk,filenami(1:ilast),indx,ztail,ier)
          if(zpath.eq.' ') then
             filename=trim(ztail)
          else
             filename=trim(zpath)//trim(ztail)
          endif
       endif
    else
       filename=filenami
    endif
    if(ier.ne.0) return
 
    mtest=filename(1:4)
    call uupper(mtest)
 
    iopen = 0                ! MDS+ tree open flag
 
    if(mtest.ne.'MDS+') then
 
! prevent memory leaks, in case where G-EQdisk object is reused.
 
       call geq_free(geq,ier)
 
!  ********* file access *************
 
       open(neqdsk, file=filename, action='read', status='old', iostat=iok)
 
       if (iok /= 0) then
          ier = 1
          go to 1001
       endif
 
       geq%ident = filename
       geq%access_method = file_access
 
       read (neqdsk,2000) (geq%case_(i),i=1,6),idum,geq%nw,geq%nh
!dbg       write(6,*) geq%nw,geq%nh
       read (neqdsk,2020) geq%rdim_m, geq%zdim_m, &
            & geq%rcentr_m, geq%rleft_m, geq%zmid_m
!dbg       write(6,*) geq%rdim_m,geq%zdim_m
       read (neqdsk,2020) geq%rmaxis_m,geq%zmaxis_m, &
            & geq%simag_Wb__Rad,geq%sibry_Wb__Rad,geq%bcentr_t
!dbg       write(6,*) geq%rmaxis_m,geq%zmaxis_m
       read (neqdsk,2020) geq%current_a,zsimag, &
            & adum,geq%rmaxis_m,adum
!dbg       write(6,*) geq%current_a
       read (neqdsk,2020) geq%zmaxis_m,adum,geq%sibry_Wb__Rad,adum,adum
!dbg       write(6,*) geq%sibry_Wb__Rad
 
       if(zsimag.ne.geq%simag_Wb__Rad) then
          write(lunmsg,*) &
               'G-eqdsk file warning: simag inconsistency, 1st value taken.'
          write(lunmsg,*) '1st simag(3rd line, 3rd word)=',geq%simag_Wb__Rad
          write(lunmsg,*) '2nd simag(4th line, 2nd word)=',zsimag
       endif
 
       allocate(geq%fpol_tm(geq%nw), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
       allocate(geq%pres_Nt__M2(geq%nw), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
       allocate(geq%ffprim_T2M2Rad__Wb(geq%nw), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
       allocate(geq%pprime_NtRad__M2Wb(geq%nw), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
       allocate(geq%qpsi(geq%nw), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
 
       allocate(geq%psirz_Wb__Rad(geq%nw, geq%nh), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
 
       read (neqdsk,2020) (geq%fpol_tm(i),i=1,geq%nw)
!dbg       write(6,*) 'fpol',(geq%fpol_tm(i),i=1,geq%nw)
       read (neqdsk,2020) (geq%pres_Nt__M2(i),i=1,geq%nw)
!dbg       write(6,*) 'pres',(geq%pres_Nt__M2(i),i=1,geq%nw)
       read (neqdsk,2020) (geq%ffprim_T2M2Rad__Wb(i),i=1,geq%nw)
       read (neqdsk,2020) (geq%pprime_NtRad__M2Wb(i),i=1,geq%nw)
       read (neqdsk,2020) ((geq%psirz_Wb__Rad(i,j),i=1,geq%nw),j=1,geq%nh)
       read (neqdsk,2020) (geq%qpsi(i),i=1,geq%nw)
!dbg       write(6,*) 'qpsi',(geq%qpsi(i),i=1,geq%nw)
       read (neqdsk,2022) geq%nbbbs,geq%limitr
 
       allocate(geq%rbbbs_m(geq%nbbbs+1), geq%zbbbs_m(geq%nbbbs+1), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
       allocate(geq%rlim_m(geq%limitr+1), geq%zlim_m(geq%limitr+1), stat=iok)
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
       geq%rbbbs_m=0
       geq%zbbbs_m=0
 
       read (neqdsk,2020) (geq%rbbbs_m(i),geq%zbbbs_m(i),i=1,geq%nbbbs)
!dbg       write(6,*) 'r,zbbbs', (geq%rbbbs_m(i),geq%zbbbs_m(i),i=1,geq%nbbbs)
       call geq_ckrzc(geq%rbbbs_m(1),geq%zbbbs_m(1),geq%nbbbs,geq%nbbbs+1)
!dbg       write(6,*) 'r,z(ck)', (geq%rbbbs_m(i),geq%zbbbs_m(i),i=1,geq%nbbbs)
 
       read (neqdsk,2020) (geq%rlim_m(i),geq%zlim_m(i),i=1,geq%limitr)
!dbg       write(6,*) 'r,zlim', (geq%rlim_m(i),geq%zlim_m(i),i=1,geq%limitr)
       call geq_ckrzc(geq%rlim_m(1),geq%zlim_m(1),geq%limitr,geq%limitr+1)
 
      ! rotation data not extracted....
 
       close(neqdsk)
 
2000   format (6a8,3i4)
2020   format (5e16.9)
2022   format (2i5)
 
    else
 
! *********** MDSplus access ***************
!
!  ...analyze the specification string in "filename" argument
!
       ilen = str_length(filename)
       icln = index(filename,':')
       if((icln.le.1).or.(icln.ge.ilen)) then
          call geq_echo('MDS+ id string (no colon):',filename)
          ier=3
          go to 1000
       endif
 
       call geq_prepars(filename(1:icln-1),opt,ireduce,ier)
       if(ier.ne.0) then
          call geq_echo('MDS+ id string (header):',filename)
          ier=3
          go to 1000
       endif
 
       call trmds_parse(filename(icln+1:ilen),mds_server,mds_tree, &
            & mds_id1,mds_id2,mds_tstr,ier)
       if((ier.ne.0).or.(mds_id2.ne.' ')) then
          call geq_echo('MDS+ id string:',filename)
          ier=3
          go to 1000
       endif
 
       if(mds_tstr.eq.' ') then
          call geq_echo('MDS+ id string (missing time):',filename)
          ier=4
          go to 1000
       endif
 
       read(mds_tstr,'(G20.0)',iostat=ier) ztime
       if(ier.ne.0) then
          call geq_echo('MDS+ id string (bad time):',filename)
          ier=4
          go to 1000
       endif
 
       ident=' '
       if(opt.eq.' ') opt='NONE'
       ilen=str_length(opt)
 
       write(ident,'(''MDS+,OPT='',a,'',REDUCE='',L1)') opt(1:ilen),ireduce
 
       ict=str_length(ident)
       ident=ident(1:ict)//':'//mds_server
 
       ict=str_length(ident)
       ident=ident(1:ict)//':'//mds_tree
 
       ict=str_length(ident)
       ilen=str_length(mds_id1)
       ident=ident(1:ict)//'('//mds_id1(1:ilen)//')'
 
       ilen=str_length(mds_id1)
       z20=' '
       z20(20-ilen+1:20)=mds_id1(1:ilen)
       read(z20,'(I20)',iostat=ier) ishot
       if(ier.ne.0) then
          call geq_echo('MDS+ id string (bad shot#):',filename)
          ier=3
          go to 1000
       endif
 
       iopen=0
 
       if(ident.ne.geq%ident) then
!  new MDSplus tree
          call geq_free(geq,ier)
          geq%ident = ident
          geq%mdsopt = opt
          geq%mdslice = ireduce
          geq%mdstime = ztime
          geq%mds_shot = ishot
          if(ireduce) then
             geq%access_method = mds_slice_access
          else
             geq%access_method = mds_full_access
          endif
          iopen=1
          call geq_mdsopen(mds_server,mds_tree,ishot,ier)
          if(ier.ne.0) then
             call geq_echo('MDS+ tree connect/open failure:',filename)
             ier=6
             go to 1000
          endif
!  read EFIT tree information
          call geq_mds_basic(geq,mds_server,mds_tree,ier)
 
!  if the full tree contains only one time point, the internal access method
!  is modified to slice access...

          ireduce = (geq%access_method .eq. mds_slice_access)

          tmplbl='xplasma:MDS+:'//geq%mdstok
          iltmp=str_length(tmplbl)
          tmplbl(iltmp+1:)=':'//mds_tree
          iltmp=str_length(tmplbl)
          iltmp2=str_length(mds_id1)
          if(iltmp+2+iltmp2.le.len(tmplbl)) then
             tmplbl(iltmp+1:)='('//mds_id1(1:iltmp2)//')'
             iltmp=str_length(tmplbl)
             iltmp2=str_length(mds_tstr)
             if(iltmp+3+iltmp2.le.len(tmplbl)) then
                tmplbl(iltmp+1:) = ' t='//mds_tstr(1:iltmp2)
             endif
          endif
 
          do i=1,6
             geq%case_(i)=tmplbl((i-1)*8+1:i*8)
          enddo
 
          if(ier.ne.0) go to 1000
       else
!  old MDSplus tree, new timeslice
          call geq_free_slice(geq)
          geq%mdstime = ztime
          if(ireduce) then
             iopen=1        ! will read timeslice from tree
             call geq_mdsopen(mds_server,mds_tree,ishot,ier)
             if(ier.ne.0) then
                call geq_echo('MDS+ tree connect/open failure:',filename)
                ier=6
                go to 1000
             endif
          endif
       endif
!
!  OK get the current timeslice
!
       call geq_mds_slice(geq,mds_tree,ier)
!
!  all done
!
    endif
 
1000 continue
 
    if(ier.eq.0) then
       geq%marker = 'xplasma:ok'
    else
       geq%marker = 'xplasma:xx'
       geq%ident = ' '
    endif
 
    if(iopen.eq.1) then
       call geq_mdsclose(mds_server,mds_tree,ishot,iok)
       iopen=0
       if(iok.ne.0) then
          call geq_echo('MDS+ tree disconnect/close failure:',filename)
          ier=6
       endif
    endif

    ! check that bdy is entirely contained inside the grid domain--
    ! flag an error if not.

    call geq_bchk(geq, 'plasma_bdy(RBBBS_m,ZBBBS_m)', &
         geq%nbbbs, geq%rbbbs_m, geq%zbbbs_m, ier)

    return
 
1001 continue
 
    ! failed to open the geqdsk file
 
    write(lunmsg,*) ' ?geq_init: file ', filename,' does not exist'
    geq%marker = 'xplasma:xx'
    geq%ident = ' '
    return
 
  end subroutine geq_init
!------------------------------------------------------
  subroutine geq_mdsopen(mds_server,mds_tree,ishot,ier)
 
! connect to server, open tree
 
    implicit NONE
 
    character*(*),intent(in) :: mds_server   ! server name or "LOCAL"
    character*(*),intent(in) :: mds_tree     ! tree name
    integer,intent(in) :: ishot              ! shot number
 
    integer,intent(out) :: ier               ! status code, 0=OK
 
    integer mds_open,mdsfconnect
 
! local variables
 
    integer lunout
    integer socket
    integer str_length,ilen,ilens
    integer istat
    character*64 mds_server_use
 
! executable code
 
    ier=0
    call geqdsk_lunmsg_get(lunout)
 
    if((mds_server.ne.'LOCAL').and.(mds_server.ne.'CONNECTED')) then
 
       write(lunout,*) ' %geqdsk_mod:  connecting to server:  ',mds_server
       mds_server_use=mds_server
       ilen=str_length(mds_server_use)+1
       mds_server_use(ilen:ilen)=char(0)
       ilen=ilen-1
 
       socket = mdsfconnect(mds_server_use)
 
       if(socket.eq.0) then
          write(lunout,*) ' ?geq_mdsopen:  failed to connect to server:', &
               mds_server_use(ilen:ilen)
          ier=1
          return
       endif
    endif
 
    if(mds_tree.ne.'OPENED') then
       write(lunout,*) ' %geqdsk_mod:  opening EFIT tree:  ',mds_tree,ishot
       istat=mds_open(mds_tree,ishot)
       if(mod(istat,2).ne.1) then
          call mdserr(lunout,'?geq_mdsopen: mds_open:  ',istat)
          ier=1
       else
          ilens=str_length(mds_server)
          ilen=str_length(mds_tree)
       endif
    endif
 
    return
  end subroutine geq_mdsopen
!------------------------------------------------------
  subroutine geq_mdsclose(mds_server,mds_tree,ishot,ier)
 
! close tree, disconnect from server
 
    implicit NONE
 
    character*(*),intent(in) :: mds_server   ! server name or "LOCAL"
    character*(*),intent(in) :: mds_tree     ! tree name
    integer,intent(in) :: ishot              ! shot number
 
    integer,intent(out) :: ier               ! status code, 0=OK
 
    integer mds_close
 
! local variables
 
    integer lunout
    integer istat
 
! executable code
 
    ier=0
    call geqdsk_lunmsg_get(lunout)
 
    if(mds_tree.ne.'OPENED') then
       istat=mds_close(mds_tree,ishot)
       if(mod(istat,2).ne.1) then
          call mdserr(lunout,'?geq_mdsclose: mds_close:  ',istat)
          ier=1
       endif
    endif
 
    if((mds_server.ne.'LOCAL').and.(mds_server.ne.'CONNECTED')) then
       call MdsCacheUnconnect  ! invalidate mdsplus socket but do not close connection
    endif
 
  end subroutine geq_mdsclose
!------------------------------------------------------
!  read EFIT MDS+ tree for all data except psi(R,Z,t)
!
!  ***dmc*** tested on NSTX only, units conversion issues & array
!            ordering issues for general EFIT trees not yet addressed.
!
  subroutine geq_mds_basic(geq,mds_server,mds_tree,ier)
 
    use geqdsk_aux
    use cgeq_paths
 
    implicit NONE
 
    type(geqdsk), intent(inout) :: geq
    character*(*), intent(in) :: mds_server
    character*(*), intent(in) :: mds_tree
    integer, intent(out) :: ier
 
! local
 
    integer ifirst,nt_orig,nt_act
    integer nh_orig,nh_min,nh_max
    integer nw_orig,nw_min,nw_max
    integer istat,istat2
    integer ilt,iret,str_length,indxs
    integer :: irank,idim3(3),insingl,inns,inn1

    integer i,j,num,ilen,iok,idum,mds_value,iwarn,idata,ick
 
    integer idims(2),idescr_long,idescr_floatarr,idescr_cstring
 
    integer, dimension(:), allocatable :: ibuf1,ibuf2
    real(geq_r8), dimension(:), allocatable :: zbuf1,ztime,zpsi0,zpsia,zcur
    real, dimension(:), allocatable :: zbuflim
    real, dimension(:), allocatable :: r4buf
    real, dimension(:,:), allocatable :: zlim2

    logical, dimension(:), allocatable :: nanflag

    real(geq_r8), dimension(:,:,:), allocatable :: zbuf3
    real(geq_r8) zero
 
    real(auxr8) factor
 
    character*80 expr,prefix,suffix,expr2
    character*1 bksl
    character*10 zmach,rlim_node,zlim_node
    character*6 ztest6
 
! executable code
 
    call ctable_init  ! initialize units conversion tables
 
    ier=0
    bksl = char(92)   ! backslash
 
    zero = 0
 
    ilt=str_length(mds_tree)
 
    write(lunmsg,*) &
         & ' %geqdsk_mod:  EFIT MDSplus tree: timebase & scalar reads.'
 
! tokamak name
 
    expr = 'machine()'
    zmach = ' '
    istat = mds_value(expr,idescr_cstring(zmach),iret)
    if(mod(istat,2).ne.1) go to 999
 
    geq%mdstok=zmach
 
    call uupper(zmach)
 
    rlim_node='xlim'
    zlim_node='ylim'
 
    if((zmach.ne.'CMOD').and.(zmach.ne.'NSTX').and. &
         (zmach.ne.'D3D')) then
       write(lunmsg,*) ' mds_value("machine()") returned "',trim(zmach),'".'
       ick=0
       if(index(mds_server,'ATLAS.GAT.COM').gt.0) then
          write(lunmsg,*) ' mds_server: "',trim(mds_server),'"; D3D assumed.'
          zmach='D3D'
          ick=1
       else if(index(mds_server,'CMODA.PSFC.MIT.EDU').gt.0) then
          write(lunmsg,*) ' mds_server: "',trim(mds_server),'"; CMOD assumed.'
          zmach='CMOD'
          ick=1
       else if(index(mds_server,'.PPPL.GOV:8501').gt.0) then
          write(lunmsg,*) ' mds_server: "',trim(mds_server),'"; NSTX assumed.'
          zmach='NSTX'
          ick=1
       else
          if(zmach.eq.' ') zmach='UNKNOWN'
       endif
       if(ick.eq.0) then
          write(lunmsg,*) ' failed to infer machine from server name: ', &
               mds_server
          write(lunmsg,*) ' default settings used for EFIT tree reads-- ', &
               'ERRORS are likely.'
       endif
    endif
 
    if(zmach.eq.'CMOD') then
       geq%kt_sw2=1
       geq%kt_sw3=3
       cgeq_apath='efit_aeqdsk:'
       cgeq_gpath='efit_geqdsk:'
       cgeq_nbdry='nbbbs'
    else if(zmach.eq.'NSTX') then
       geq%kt_sw2=2
       geq%kt_sw3=3
       ztest6 = mds_tree(1:6)
       call uupper(ztest6)
       if(ztest6.eq.'LRDFIT') then
          cgeq_apath='lrdfit_aeqdsk:'
          cgeq_gpath='lrdfit_geqdsk:'
       else
          cgeq_apath='efit_aeqdsk:'
          cgeq_gpath='efit_geqdsk:'
       endif
       cgeq_nbdry='nbbbs'
    else if(zmach.eq.'D3D') then
       geq%kt_sw2=2
       geq%kt_sw3=3
       cgeq_apath=' '
       cgeq_gpath=' '
       cgeq_nbdry='nbdry'
    else
       write(lunmsg,*) ' **geqdsk_mod warning:  EFIT tree for machine ',zmach
       write(lunmsg,*) '   unknown machine, profile storage ordering guessed.'
       geq%kt_sw2=1
       geq%kt_sw3=3
    endif
 
! #of time pts
 
    call geq_mk_expr('size(','gtime',')',expr)
    istat = mds_value(expr,idescr_long(nt_orig),iret)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(r4buf(nt_orig),stat=ier)
    if(ier.ne.0) go to 995
    allocate(ibuf1(nt_orig),stat=ier)
    if(ier.ne.0) go to 995
    allocate(zbuf1(nt_orig),stat=ier)
    if(ier.ne.0) go to 995
    allocate(ztime(nt_orig),zpsi0(nt_orig),zpsia(nt_orig),stat=ier)
    if(ier.ne.0) go to 995
    allocate(nanflag(nt_orig),stat=ier); nanflag = .FALSE.
    if(ier.ne.0) go to 995
    allocate(zcur(nt_orig),stat=ier)
    if(ier.ne.0) go to 995
 
    idims(1)=nt_orig
    idims(2)=0
 
! EFIT tree timebase
 
    call geq_mk_expr('data(','gtime',')',expr)
    istat = mds_value(expr,idescr_floatarr(r4buf(1),idims,1),iret)
    if(mod(istat,2).ne.1) go to 999
    call geq_nanscan(expr,nt_orig,r4buf,nanflag)
    ztime=r4buf
 
    call geq_units_check(expr,ztime,nt_orig,time_units)
 
! psi(0), psi(a), plasma current
 
    call geq_mk_expr('data(','ssimag',')',expr)
    istat = mds_value(expr,idescr_floatarr(r4buf(1),idims,1),iret)
    if(mod(istat,2).ne.1) go to 999
    call geq_nanscan(expr,nt_orig,r4buf,nanflag)
    zpsi0=r4buf
 
    call geq_units_check(expr,zpsi0,nt_orig,psi_units)
 
    call geq_mk_expr('data(','ssibry',')',expr)
    istat = mds_value(expr,idescr_floatarr(r4buf(1),idims,1),iret)
    if(mod(istat,2).ne.1) go to 999
    call geq_nanscan(expr,nt_orig,r4buf,nanflag)
    zpsia=r4buf
 
    call geq_units_check(expr,zpsia,nt_orig,psi_units)
 
    call geq_mk_expr('data(','cpasma',')',expr)
    istat = mds_value(expr,idescr_floatarr(r4buf(1),idims,1),iret)
    if(mod(istat,2).ne.1) go to 999
    call geq_nanscan(expr,nt_orig,r4buf,nanflag)
    zcur=r4buf
 
    call geq_units_check(expr,zcur,nt_orig,current_units)
 
! DMC -- for NSTX, add this read for NaNscan purposes only...
    call geq_mk_expr('data(','WMHD',')',expr)
    indxs=index(expr,'_geqdsk')+1
    if(indxs.le.1) then
       write(lunmsg,*) ' %note-- no path to _aeqdsk; no WMHD read attempt.'
       r4buf=0.0
    else
       expr(indxs:indxs)='a'  ! _aeqdsk subtree...
       istat = mds_value(expr,idescr_floatarr(r4buf(1),idims,1),iret)
       if(mod(istat,2).ne.1) then
          write(lunmsg,*) ' %warning, MDSplus read failed on: ',trim(expr), &
               ' (OK to ignore).'
          r4buf=0.0
       endif
       call geq_nanscan(expr,nt_orig,r4buf,nanflag)
    endif
 
! assemble list of valid times:  cur.ne.0, psi0.ne.psia
 
    geq%curtime = 0

    allocate(geq%gtime_s(nt_orig),stat=ier)
    if(ier.ne.0) go to 995
    allocate(geq%t_index(nt_orig),stat=ier)
    if(ier.ne.0) go to 995
 
    nt_act = 0
 
    ifirst=0
    do i=1,nt_orig
       geq%gtime_s(i)=zero
       if((.not.nanflag(i)).and. &
            (zcur(i).ne.zero).and.((zpsia(i)-zpsi0(i)).ne.zero)) then
          if(ifirst.eq.0) ifirst=i
          nt_act = nt_act + 1
          geq%t_index(i) = nt_act
          geq%gtime_s(nt_act) = ztime(i)
       else
          geq%t_index(i)=0
       endif
    enddo
 
    if(ifirst.eq.0) then
       write(lunmsg,*) ' ?geqdsk_mod:  no valid time points in EFIT tree.'
       ier=7
       return
    endif
 
    write(lunmsg,'(a,2(1pe12.5,1x))') &
         & ' %geqdsk:  available time range (s): ', &
         & geq%gtime_s(1),geq%gtime_s(nt_act)
 
! save current, psi0, psia in "geq".
 
    geq%t_num_valid = nt_act
    if(nt_act.le.1) then
       geq%access_method = mds_slice_access
    endif

    allocate(geq%t_current(nt_act),geq%t_simag(nt_act),geq%t_sibry(nt_act), &
         & stat=ier)
    if(ier.ne.0) go to 995
 
    do i=1,nt_orig
       j=geq%t_index(i)
       if(j.ne.0) then
          geq%t_current(j)=zcur(i)
          geq%t_simag(j)=zpsi0(i)
          geq%t_sibry(j)=zpsia(i)
       endif
    enddo
 
! read the remaining scalar fcns of time, save in "geq"
 
    allocate(geq%t_rdim(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_rdim(1),nt_act, &
         & geq%t_index(1),'xdim', &
         & mds_tree,expr,length_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
!
! check for no time variation in xdim and related quantities:
!
    geq%kt_max=nt_act
    if(geq%t_rdim(2).le.0) geq%kt_max=1
 
    allocate(geq%t_zdim(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_zdim(1),nt_act, &
         geq%t_index(1),'zdim', &
         & mds_tree,expr,length_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(geq%t_rcentr(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_rcentr(1),nt_act, &
         & geq%t_index(1),'rbcent', &
         & mds_tree,expr,length_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(geq%t_rleft(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_rleft(1),nt_act, &
         & geq%t_index(1),'rgrid1', &
         & mds_tree,expr,length_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(geq%t_zmid(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_zmid(1),nt_act, &
         & geq%t_index(1),'zmid', &
         & mds_tree,expr,length_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(geq%t_rmaxis(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_rmaxis(1),nt_act, &
         & geq%t_index(1),'rmaxis', &
         & mds_tree,expr,length_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(geq%t_zmaxis(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_zmaxis(1),nt_act, &
         & geq%t_index(1),'zmaxis', &
         & mds_tree,expr,length_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(geq%t_bcentr(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_sub1(zbuf1,nt_orig,geq%t_bcentr(1),nt_act, &
         geq%t_index(1),'bcentr', &
         & mds_tree,expr,bfield_units,istat)
    if(mod(istat,2).ne.1) go to 999
 
! get the spatial grid dimensions; require that these not change as a
! function of time
 
    allocate(geq%t_nw(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_subi(ibuf1,nt_orig,geq%t_nw(1),nt_act,geq%t_index(1),'mw', &
         & mds_tree,expr,istat)
    if(mod(istat,2).ne.1) go to 999
 
    allocate(geq%t_nh(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_subi(ibuf1,nt_orig,geq%t_nh(1),nt_act,geq%t_index(1),'mh', &
         & mds_tree,expr,istat)
    if(mod(istat,2).ne.1) go to 999
 
    nw_orig=geq%t_nw(1)
    nh_orig=geq%t_nh(1)
 
    nw_min=nw_orig
    nw_max=nw_orig
    nh_min=nh_orig
    nh_max=nh_orig
 
! check that there are no changes in time.
 
    ier = 0
    do i=1,geq%kt_max
       nh_min=min(nh_min,geq%t_nh(i))
       nh_max=max(nh_max,geq%t_nh(i))
       nw_min=min(nw_min,geq%t_nw(i))
       nw_max=max(nw_max,geq%t_nw(i))
    enddo
    if(nh_min.ne.nh_max) then
       ier=7
       write(lunmsg,*) ' ?geq_mds_basic:  gridsize MH not constant in time.'
       write(lunmsg,*) '  min(MH) = ',nh_min,' max(MH) = ',nh_max
    endif
    if(nw_min.ne.nw_max) then
       ier=7
       write(lunmsg,*) ' ?geq_mds_basic:  gridsize MW not constant in time.'
       write(lunmsg,*) '  min(MW) = ',nw_min,' max(MW) = ',nw_max
    endif
 
! set some constants in the main (non-time-varying) part of geq
 
    geq%nw = nw_orig
    geq%nh = nh_orig
 
! read limiter & boundary data
    write(lunmsg,*) &
         & ' %geqdsk_mod:  EFIT MDSplus tree: limiter & bdy reads.'
 
    allocate(ibuf2(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_subi(ibuf1,nt_orig,ibuf2,nt_act,geq%t_index(1),'limitr', &
         & mds_tree,expr,istat)
    if(mod(istat,2).ne.1) go to 999
 
    geq%limitr = ibuf2(1)   ! no. of limiter pts, assumed no change vs. time
 
    allocate(geq%rlim_m(geq%limitr+1),geq%zlim_m(geq%limitr+1),stat=ier)
    if(ier.ne.0) go to 995
    allocate(zbuflim(geq%limitr),stat=ier)
    if(ier.ne.0) go to 995
 
    geq%rlim_m = 0
    geq%zlim_m = 0
 
! read limiter R pts
 
    iwarn=1
 
    call geq_mk_expr('rank(','lim',')',expr)
    istat = mds_value(expr,idescr_long(idata),iret)
 
    if((mod(istat,2).eq.1).and.(idata.eq.2)) then
       call geq_mk_expr('size(','lim',')',expr)
       istat = mds_value(expr,idescr_long(idata),iret)
 
       if(mod(istat,2).eq.1) then
          if(idata.ne.2*geq%limitr) then
             write(lunmsg,*) &
                  ' %geqdsk_mod:  size(lim).ne.2*limitr; attempting recovery.'
          else
             iwarn=0
          endif
       endif
    endif
 
    if(iwarn.eq.0) then
 
! read "lim" node and unpack into rlim_m, zlim_m arrays -- no time variation
 
       allocate(zlim2(2,geq%limitr),stat=ier)
       if(ier.ne.0) go to 995
 
       idims(1)=2
       idims(2)=geq%limitr
       call geq_mk_expr('data(','lim',')',expr)
       istat = mds_value(expr,idescr_floatarr(zlim2(1,1),idims,2),iret)
       if(mod(istat,2).ne.1) go to 999
 
       geq%rlim_m(1:geq%limitr)=zlim2(1,1:geq%limitr)
       geq%zlim_m(1:geq%limitr)=zlim2(2,1:geq%limitr)
 
       call geq_units_check(expr,geq%rlim_m(1),geq%limitr,length_units)
       call geq_units_check(expr,geq%zlim_m(1),geq%limitr,length_units)
 
    else
 
! read separate nodes for R and Z; may or may not seem to have time variation
 
       idims(1)=geq%limitr
       idims(2)=0
 
! read limiter R pts

!   DMC -- code that follows deals with some xlim // ylim LRDFIT data
!          that are nominally 3d arrays
       insingl=0
       inns=0
       inn1=0

       call geq_mk_expr('rank(',rlim_node,')',expr)
       istat = mds_value(expr,idescr_long(irank),iret)

       if(irank.eq.3) then
          call geq_mk_expr('size(',rlim_node,',0)',expr)
          istat = mds_value(expr,idescr_long(idata),iret)
          idim3(1)=idata
          if(idata.eq.1) then
             insingl=insingl+1
             inn1=1
          else
             inns=1
          endif

          call geq_mk_expr('size(',rlim_node,',1)',expr)
          istat = mds_value(expr,idescr_long(idata),iret)
          idim3(2)=idata
          if(idata.eq.1) then
             insingl=insingl+1
             inn1=2
          else
             inns=2
          endif

          call geq_mk_expr('size(',rlim_node,',2)',expr)
          istat = mds_value(expr,idescr_long(idata),iret)
          idim3(3)=idata
          if(idata.eq.1) then
             insingl=insingl+1
             inn1=3
          else
             inns=3
          endif

          if(insingl.eq.2) then
             if(inns.eq.1) then
                suffix = ')[*,0,0]'
             else if(inns.eq.2) then
                suffix = ')[0,*,0]'
             else if(inns.eq.3) then
                suffix = ')[0,0,*]'
             endif
          else if(insingl.eq.1) then
             if(inn1.eq.1) then
                if(geq%kt_sw2.eq.2) then
                   write(suffix,'(")[0,*,",i6,"]")') ifirst
                else
                   write(suffix,'(")[0,",i6,",*]")') ifirst
                endif
             else if(inn1.eq.2) then
                if(geq%kt_sw2.eq.2) then
                   write(suffix,'(")[*,0,",i6,"]")') ifirst
                else
                   write(suffix,'(")[",i6,",0,*]")') ifirst
                endif
             else
                if(geq%kt_sw2.eq.2) then
                   write(suffix,'(")[*,",i6,",0]")') ifirst
                else
                   write(suffix,'(")[",i6,",*,0]")') ifirst
                endif
             endif
          endif
       endif

       if((insingl.eq.0).or.(insingl.eq.3)) then
          ! normal rank=2 xlim & ylim nodes pass through here...

          if(geq%kt_sw2.eq.2) then
             write(suffix,'(")[*,",i6,"]")') ifirst
          else
             write(suffix,'(")[",i6,",*]")') ifirst
          endif
       endif

       call geq_mk_expr('data(',rlim_node,suffix,expr)
 
       istat = mds_value(expr,idescr_floatarr(zbuflim(1),idims,1),iret)
       if(mod(istat,2).ne.1) then
          expr2=expr
          istat2=istat
          suffix = ')[*]'
! try for 1d object (no time dependence)
          call geq_mk_expr('data(',rlim_node,suffix,expr)
          istat = mds_value(expr,idescr_floatarr(zbuflim(1),idims,1),iret)
          if(mod(istat,2).ne.1) then
             call mdserr(lunmsg,'?geq_mds_basic: tried: '//expr2,istat2)
             go to 999
          endif
       endif
       geq%rlim_m(1:geq%limitr)=zbuflim
       call geq_units_check(expr,geq%rlim_m(1),geq%limitr,length_units)
 
! read limiter Z pts
 
       call geq_mk_expr('data(',zlim_node,suffix,expr)
 
       istat = mds_value(expr,idescr_floatarr(zbuflim(1),idims,1),iret)
       if(mod(istat,2).ne.1) then
          expr2=expr
          istat2=istat
          suffix = ')[*]'
! try for 1d object (no time dependence)
          call geq_mk_expr('data(',zlim_node,suffix,expr)
          istat = mds_value(expr,idescr_floatarr(zbuflim(1),idims,1),iret)
          if(mod(istat,2).ne.1) then
             call mdserr(lunmsg,'?geq_mds_basic: tried: '//expr2,istat2)
             go to 999
          endif
       endif
       geq%zlim_m(1:geq%limitr)=zbuflim
       call geq_units_check(expr,geq%zlim_m(1),geq%limitr,length_units)
 
    endif
 
! make sure points are separate and that the limiter contour closes
    call geq_ckrzc(geq%rlim_m(1),geq%zlim_m(1),geq%limitr,geq%limitr+1)
 
! read plasma bdy data
 
! no. of bdy pts, vs. time
    allocate(geq%t_nbbbs(nt_act),stat=ier)
    if(ier.ne.0) go to 995
    call geq_mds_subi(ibuf1,nt_orig,geq%t_nbbbs(1),nt_act, &
         & geq%t_index(1),cgeq_nbdry, &
         & mds_tree,expr,istat)
    if(mod(istat,2).ne.1) go to 999
 
    call geq_mk_expr('rank(','bdry',')',expr)
    istat = mds_value(expr,idescr_long(idata),iret)
    if(mod(istat,2).ne.1) idata=0
!
! now deciding whether to read all time dependent information...
!
    if(geq%access_method.eq.mds_full_access) then
 
! This block only reached if there are 2 or more valid time points

! array data for:
!   R plasma bdy pts
!   Z plasma bdy pts
 
       num=maxval(geq%t_nbbbs)
 
       allocate(geq%t_rbbbs(num,nt_act),stat=ier)
       if(ier.ne.0) go to 995
 
       allocate(geq%t_zbbbs(num,nt_act),stat=ier)
       if(ier.ne.0) go to 995
 
       if((geq%kt_sw2.eq.2).and.(idata.eq.3)) then
! get the data, combined, from the "bdry" node.
 
          allocate(zbuf3(2,num,nt_act),stat=ier)
          if(ier.ne.0) go to 995
 
          call geq_mds_sub3(2,num,nt_orig,zbuf3,nt_act,3, &
               & geq%t_index(1),'bdry', &
               & mds_tree,expr,length_units,istat)
          if(mod(istat,2).ne.1) go to 999
 
          geq%t_rbbbs = zbuf3(1,1:num,1:nt_act)
          geq%t_zbbbs = zbuf3(2,1:num,1:nt_act)
 
          deallocate(zbuf3)
 
! get the data, in a conventional way, from rbbbs & zbbbs
 
       else
          call geq_mds_sub2(num,nt_orig,geq%t_rbbbs(1,1),nt_act, geq%kt_sw2, &
               & geq%t_index(1),'rbbbs', &
               & mds_tree,expr,length_units,istat)
          if(mod(istat,2).ne.1) go to 999
 
          call geq_mds_sub2(num,nt_orig,geq%t_zbbbs(1,1),nt_act, geq%kt_sw2, &
               & geq%t_index(1),'zbbbs', &
               & mds_tree,expr,length_units,istat)
          if(mod(istat,2).ne.1) go to 999
       endif
 
       write(lunmsg,*) &
            & ' %geqdsk_mod:  EFIT MDSplus tree: f(psi,t) profile reads.'
 
       num=geq%nw
 
! f
       allocate(geq%t_fpol(num,nt_act),stat=ier)
       if(ier.ne.0) go to 995
       call geq_mds_sub2(num,nt_orig,geq%t_fpol(1,1),nt_act, geq%kt_sw2, &
            & geq%t_index(1),'fpol', &
            & mds_tree,expr,fpol_units,istat)
       if(mod(istat,2).ne.1) go to 999
 
! P
       allocate(geq%t_pres(num,nt_act),stat=ier)
       if(ier.ne.0) go to 995
       call geq_mds_sub2(num,nt_orig,geq%t_pres(1,1),nt_act, geq%kt_sw2, &
            & geq%t_index(1),'pres', &
            & mds_tree,expr,pressure_units,istat)
       if(mod(istat,2).ne.1) go to 999
 
! ffprime
       allocate(geq%t_ffprim(num,nt_act),stat=ier)
       if(ier.ne.0) go to 995
       call geq_mds_sub2(num,nt_orig,geq%t_ffprim(1,1),nt_act, geq%kt_sw2, &
            & geq%t_index(1),'ffprim',&
            & mds_tree,expr,ffprime_units,istat)
       if(mod(istat,2).ne.1) go to 999
 
! PPrime
       allocate(geq%t_pprime(num,nt_act),stat=ier)
       if(ier.ne.0) go to 995
       call geq_mds_sub2(num,nt_orig,geq%t_pprime(1,1),nt_act, geq%kt_sw2, &
            & geq%t_index(1),'pprime',&
            & mds_tree,expr,pprime_units,istat)
       if(mod(istat,2).ne.1) go to 999
 
! q
       allocate(geq%t_qpsi(num,nt_act),stat=ier)
       if(ier.ne.0) go to 995
       call geq_mds_sub2(num,nt_orig,geq%t_qpsi(1,1),nt_act, geq%kt_sw2, &
            & geq%t_index(1),'qpsi', &
            & mds_tree,expr,no_units,istat)
       if(mod(istat,2).ne.1) go to 999
 
 
       write(lunmsg,*) &
            & ' %geqdsk_mod:  EFIT MDSplus tree: psi(R,Z,t) profile read.'
 
       allocate(geq%t_psirz(geq%nw,geq%nh,nt_act),stat=ier)
       if(ier.ne.0) go to 995
       call geq_mds_sub3(geq%nw,geq%nh,nt_orig,geq%t_psirz(1,1,1),nt_act, &
            & geq%kt_sw3, geq%t_index(1),'psirz', &
            & mds_tree,expr,psi_units,istat)
       if(mod(istat,2).ne.1) go to 999
 
       write(lunmsg,*) ' %geqdsk_mod:  read completed.'
 
    endif
 
    go to 1000
 
! ************** errors **************
 
995 continue
    write(lunmsg,*) ' ?allocation error in geq_mds_basic:  ',ier
    ier = 2
    go to 1000
 
999 continue
    call mdserr(lunmsg,'?geq_mds_basic: i/o: '//expr,istat)
    ier = 6
    go to 1000
 
1000 continue
 
    deallocate(ibuf1,r4buf,zbuf1,ztime,zpsi0,zpsia,zcur,stat=idum)
    deallocate(ibuf2,stat=idum)
    deallocate(zbuflim,stat=idum)
    deallocate(zlim2,stat=idum)
    deallocate(nanflag,stat=idum)
 
    return
 
  end subroutine geq_mds_basic
!------------------------------------------------------
  subroutine geq_mds_slice(geq,mds_tree,ier)
 
    use geqdsk_aux
 
    implicit NONE
 
    type(geqdsk), intent(inout) :: geq
    character*(*), intent(in) :: mds_tree
    integer, intent(out) :: ier
 
! local
 
    integer i,imin,imin1,imin2,iminf,iok,istat,iret,idata
    integer inb1,inb2
    real(geq_r8) dtmin,dt,zztime,zf1,zf2
 
    integer :: inbbb

    integer ilt,str_length
    integer idims(2)
    character*80 expr,prefix,suffix,suffix2
    character*1 bksl
 
    integer mds_value,idescr_floatarr,idescr_long
 
    real, dimension(:), allocatable :: zbuf1,zbuf1a
    real, dimension(:,:), allocatable :: zbuf2
 
! executable code
! ... find nearest timepoint
 
    imin=1
    dtmin=abs(geq%mdstime-geq%gtime_s(1))
    do i=1,geq%t_num_valid
       dt=abs(geq%mdstime-geq%gtime_s(i))
       if(dt.lt.dtmin) then
          imin=i
          dtmin=dt
       endif
    enddo

    if((geq%access_method.eq.mds_full_access).and. &
         (geq%mdstime.gt.geq%gtime_s(1)).and. &
         (geq%mdstime.lt.geq%gtime_s(geq%t_num_valid))) then
       ! time interpolation OK
       if(geq%mdstime.ge.geq%gtime_s(imin)) then
          imin1=imin
          imin2=imin+1
       else
          imin1=imin-1
          imin2=imin
       endif
       zf2=(geq%mdstime-geq%gtime_s(imin1))/ &
            (geq%gtime_s(imin2)-geq%gtime_s(imin1))
       zf1=1.0_geq_r8-zf2
       write(lunmsg, &
        '('' %geqdsk_mod:  '',1pe12.5,'' btw 2 EFIT times: '',2(1pe12.5,1x))'&
        ) &
             geq%mdstime,geq%gtime_s(imin1),geq%gtime_s(imin2)
       write(lunmsg,*) '  (times in seconds) interpolation will be done.'
       geq%curtime = geq%mdstime
       geq%nbbbs = min(401,(4*maxval(geq%t_nbbbs(1:geq%t_num_valid))))
    else
       ! no time interpolation: either not full access, or, time out
       ! of range...
       imin1=imin
       imin2=imin
       zf1=1.0_geq_r8
       zf2=0.0_geq_r8
       write(lunmsg,'('' %geqdsk_mod:  selected EFIT time:  '',1pe12.5)') &
            & geq%gtime_s(imin)
       geq%curtime = geq%gtime_s(imin)
       geq%nbbbs = geq%t_nbbbs(imin)
    endif

    geq%rdim_m = zf1*geq%t_rdim(min(imin1,geq%kt_max)) + &
         zf2*geq%t_rdim(min(imin2,geq%kt_max))

    geq%zdim_m = zf1*geq%t_zdim(min(imin1,geq%kt_max)) + &
         zf2*geq%t_zdim(min(imin2,geq%kt_max))

    geq%rcentr_m = zf1*geq%t_rcentr(min(imin1,geq%kt_max)) + &
         zf2*geq%t_rcentr(min(imin2,geq%kt_max))

    geq%rleft_m = zf1*geq%t_rleft(min(imin1,geq%kt_max)) + &
         zf2*geq%t_rleft(min(imin2,geq%kt_max))

    geq%zmid_m = zf1*geq%t_zmid(min(imin1,geq%kt_max)) + &
         zf2*geq%t_zmid(min(imin2,geq%kt_max))

    geq%rmaxis_m = zf1*geq%t_rmaxis(imin1) + &
         zf2*geq%t_rmaxis(imin2)
    geq%zmaxis_m = zf1*geq%t_zmaxis(imin1) + &
         zf2*geq%t_zmaxis(imin2)
    geq%bcentr_t = zf1*geq%t_bcentr(imin1) + &
         zf2*geq%t_bcentr(imin2)
    geq%current_a = zf1*geq%t_current(imin1) + &
         zf2*geq%t_current(imin2)
    geq%simag_Wb__Rad = zf1*geq%t_simag(imin1) + &
         zf2*geq%t_simag(imin2)
    geq%sibry_Wb__Rad = zf1*geq%t_sibry(imin1) + &
         zf2*geq%t_sibry(imin2)
    
    allocate(zbuf1(geq%nw), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
    allocate(geq%fpol_tm(geq%nw), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
    allocate(geq%pres_Nt__M2(geq%nw), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
    allocate(geq%ffprim_T2M2Rad__Wb(geq%nw), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
    allocate(geq%pprime_NtRad__M2Wb(geq%nw), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
    allocate(geq%qpsi(geq%nw), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
 
    allocate(geq%psirz_Wb__Rad(geq%nw, geq%nh), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
 
! here, allow extra word in case endpt connection needs to be made
 
    allocate(zbuf1a(geq%nbbbs+1), stat=iok); zbuf1a=0.0
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
    allocate(geq%rbbbs_m(geq%nbbbs+1),geq%zbbbs_m(geq%nbbbs+1), stat=iok)
    if(iok /= 0) then
       ier = 2
       go to 1000
    endif
 
    geq%rbbbs_m=0
    geq%zbbbs_m=0
 
! OK, either copy the slices from data already read, or, read the data now.
 
    if(geq%access_method.eq.mds_full_access) then
 
       geq%fpol_tm = zf1*geq%t_fpol(1:geq%nw,imin1) + &
            zf2*geq%t_fpol(1:geq%nw,imin2)
       geq%pres_Nt__M2 = zf1*geq%t_pres(1:geq%nw,imin1) + &
            zf2*geq%t_pres(1:geq%nw,imin2)
       geq%ffprim_T2M2Rad__Wb = zf1*geq%t_ffprim(1:geq%nw,imin1) + &
            zf2*geq%t_ffprim(1:geq%nw,imin2)
       geq%pprime_NtRad__M2Wb = zf1*geq%t_pprime(1:geq%nw,imin1) + &
            zf2*geq%t_pprime(1:geq%nw,imin2)
       geq%qpsi = zf1*geq%t_qpsi(1:geq%nw,imin1) + &
            zf2*geq%t_qpsi(1:geq%nw,imin2)
 
       geq%psirz_Wb__Rad = zf1*geq%t_psirz(1:geq%nw,1:geq%nh,imin1) + &
            zf2*geq%t_psirz(1:geq%nw,1:geq%nh,imin2)
 
       inb1=geq%t_nbbbs(imin1)
       inb2=geq%t_nbbbs(imin2)
       call geq_mergbb(geq%rbbbs_m(1:geq%nbbbs),geq%zbbbs_m(1:geq%nbbbs), &
            zf1,inb1,geq%t_rbbbs(1:inb1,imin1),geq%t_zbbbs(1:inb1,imin1), &
            zf2,inb2,geq%t_rbbbs(1:inb2,imin2),geq%t_zbbbs(1:inb2,imin2))
 
       ier=0
 
    else
 
! see if "bdry" node exists -- may need to use it...
 
       call geq_mk_expr('rank(','bdry',')',expr)
       istat = mds_value(expr,idescr_long(idata),iret)
       if(mod(istat,2).ne.1) idata=0
 
! read slices from tree -- presumed already open
! *** time dimension is presumed to be the last dimension ***
 
       bksl = char(92)   ! backslash
 
       ilt=str_length(mds_tree)
 
       idims(1)=1
       idims(2)=1
 
       iminf = imin                ! rediscover EFIT tree time index
       do while (geq%t_index(iminf).ne.imin)
          iminf = iminf+1
          if(iminf.gt.size(geq%gtime_s)) then
             call errmsg_exit(' ?geqdsk_mod -- internal algorithm error.')
          endif
       enddo
       iminf=iminf-1    ! *** correct to mdsplus/C-style 0 based indexing!
 
       idims(1)=geq%nw
       idims(2)=1
 
       write(lunmsg,*) &
            & ' %geqdsk_mod:  EFIT MDSplus tree: f(psi) profile read.'
 
       suffix2 = ')[*]'
       if(geq%kt_sw2.eq.2) then
          write(suffix,1001) iminf
       else
          write(suffix,1011) iminf
       endif
 
       call try_expr1('fpol')
       if(mod(istat,2).ne.1) go to 999
       geq%fpol_tm=zbuf1
       call geq_units_check(expr,geq%fpol_tm(1),geq%nw,fpol_units)
 
       call try_expr1('pres')
       if(mod(istat,2).ne.1) go to 999
       geq%pres_Nt__M2=zbuf1
       call geq_units_check(expr,geq%pres_Nt__M2(1),geq%nw,pressure_units)
 
       call try_expr1('ffprim')
       if(mod(istat,2).ne.1) go to 999
       geq%ffprim_T2M2Rad__Wb=zbuf1
       call geq_units_check(expr,geq%ffprim_T2M2Rad__Wb(1), &
            & geq%nw,ffprime_units)
 
       call try_expr1('pprime')
       if(mod(istat,2).ne.1) go to 999
       geq%pprime_NtRad__M2Wb=zbuf1
       call geq_units_check(expr,geq%pprime_NtRad__M2Wb(1),geq%nw,pprime_units)
 
       call try_expr1('qpsi')
       if(mod(istat,2).ne.1) go to 999
       geq%qpsi=zbuf1
       call geq_units_check(expr,geq%qpsi(1),geq%nw,no_units)
 
       if((idata.gt.0).and.(geq%kt_sw2.eq.2)) then
!  use bdry node...
          idims(1)=2
          idims(2)=geq%nbbbs
 
          write(suffix,1002) iminf
          suffix2=')[*,*]'

          allocate(zbuf2(2,geq%nbbbs), stat=iok)
 
          if(iok /= 0) then
             ier = 2
             go to 1000
          endif
 
          call try_expr2('bdry')
          if(mod(istat,2).ne.1) go to 999
 
          geq%rbbbs_m(1:geq%nbbbs) = zbuf2(1,1:geq%nbbbs)
          geq%zbbbs_m(1:geq%nbbbs) = zbuf2(2,1:geq%nbbbs)
 
          call geq_units_check(expr,geq%rbbbs_m(1),geq%nbbbs,length_units)
          call geq_units_check(expr,geq%zbbbs_m(1),geq%nbbbs,length_units)
 
          deallocate(zbuf2, stat=iok)
 
       else
!  use standard rbbbs & zbbbs nodes...
 
          idims(1)=geq%nbbbs
 
          call try_expr1a('rbbbs')
          if(mod(istat,2).ne.1) go to 999
          geq%rbbbs_m=zbuf1a
          call geq_units_check(expr,geq%rbbbs_m(1),geq%nbbbs,length_units)
 
          call try_expr1a('zbbbs')
          if(mod(istat,2).ne.1) go to 999
          geq%zbbbs_m=zbuf1a
          call geq_units_check(expr,geq%zbbbs_m(1),geq%nbbbs,length_units)
       endif
 
       call geq_ckrzc(geq%rbbbs_m(1),geq%zbbbs_m(1),geq%nbbbs,geq%nbbbs+1)
 
       write(lunmsg,*) &
            & ' %geqdsk_mod:  EFIT MDSplus tree: psi(R,Z) profile read.'
 
       idims(1)=geq%nw
       idims(2)=geq%nh
 
       allocate(zbuf2(geq%nw,geq%nh), stat=iok)
 
       if(iok /= 0) then
          ier = 2
          go to 1000
       endif
 
       suffix2=')[*,*]'
       if(geq%kt_sw3.eq.3) then
          write(suffix,1002) iminf
       else if(geq%kt_sw3.eq.1) then
          write(suffix,1012) iminf
       else
          write(suffix,1022) iminf
       endif
 
       call try_expr2('psirz')
       if(mod(istat,2).ne.1) go to 999
       geq%psirz_Wb__Rad=zbuf2
       call geq_units_check(expr,geq%psirz_Wb__Rad(1,1), &
            & geq%nw*geq%nh,psi_units)
 
       deallocate(zbuf2, stat=iok)
       ier=0
 
    endif
 
    go to 1000
 
999 continue
    call mdserr(lunmsg,'?geq_mds_slice: i/o: '//expr,istat)
    ier = 6
    go to 1000
 
1000 continue
 
1001 format(')[*,',i6,']')    ! time index last
1011 format(')[',i6,',*]')    ! time index first
1002 format(')[*,*,',i6,']')  ! time index last
1012 format(')[',i6,',*,*]')  ! time index first
1022 format(')[*,',i6,',*]')
1003 format(')[',i6,']')
 
    deallocate(zbuf1,zbuf1a,stat=iok)
    deallocate(zbuf2,stat=iok)
 
    return
 
  CONTAINS
    subroutine try_expr1(arg)
      character*(*) :: arg
      integer :: idescr_floatarr

      call geq_mk_expr('data(',arg,suffix,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf1(1),idims,1),iret)
      if(mod(istat,2).eq.1) return
      if(geq%t_num_valid.gt.1) return

      ! single time point tree -- attempt fallback

      call geq_mk_expr('data(',arg,suffix2,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf1(1),idims,1),iret)
      if(mod(istat,2).eq.1) return

      ! no -- still bad -- recover original error

      call geq_mk_expr('data(',arg,suffix,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf1(1),idims,1),iret)
    end subroutine try_expr1

    subroutine try_expr1a(arg)
      character*(*) :: arg
      integer :: idescr_floatarr

      call geq_mk_expr('data(',arg,suffix,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf1a(1),idims,1),iret)
      if(mod(istat,2).eq.1) return
      if(geq%t_num_valid.gt.1) return

      ! single time point tree -- attempt fallback

      call geq_mk_expr('data(',arg,suffix2,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf1a(1),idims,1),iret)
      if(mod(istat,2).eq.1) return

      ! no -- still bad -- recover original error

      call geq_mk_expr('data(',arg,suffix,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf1a(1),idims,1),iret)
    end subroutine try_expr1a

    subroutine try_expr2(arg)
      character*(*) :: arg
      integer :: idescr_floatarr

      call geq_mk_expr('data(',arg,suffix,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf2(1,1),idims,2),iret)
      if(mod(istat,2).eq.1) return
      if(geq%t_num_valid.gt.1) return

      ! single time point tree -- attempt fallback

      call geq_mk_expr('data(',arg,suffix2,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf2(1,1),idims,2),iret)
      if(mod(istat,2).eq.1) return

      ! no -- still bad -- recover original error

      call geq_mk_expr('data(',arg,suffix,expr)
      istat = mds_value(expr, &
           idescr_floatarr(zbuf2(1,1),idims,2),iret)
    end subroutine try_expr2

  end subroutine geq_mds_slice
 
  subroutine geq_mergbb(rbout,zbout, &
       zf1,inb1,rb1,zb1, zf2,inb2,rb2,zb2)

    real(geq_r8), dimension(:) :: rbout,zbout  ! output contour

    real(geq_r8), intent(in) :: zf1   ! weighting for 1st input contour
    integer, intent(in) :: inb1       ! 1st input contour array dimension
    real(geq_r8), intent(in) :: rb1(inb1),zb1(inb1)  ! 1st contour

    real(geq_r8), intent(in) :: zf2   ! weighting for 2nd input contour
    integer, intent(in) :: inb2       ! 2nd input contour array dimension
    real(geq_r8), intent(in) :: rb2(inb2),zb2(inb2)  ! 2nd contour

    !--------------------------------
    ! local:

    integer :: i,n1,n2,n1a,n2a,nbout
    real(geq_r8), dimension(:), allocatable :: r1,z1,r2,z2
    real(geq_r8), dimension(:), allocatable :: a1,d1,a2,d2
    real(geq_r8) :: rcen1,zcen1,rcen2,zcen2,rcen,zcen
    real(geq_r8) :: rmaxwk,rminwk
    real(geq_r8) :: zdamin,zalim,ztha,dfind1,dfind2,dfind
    !--------------------------------

    nbout=size(rbout)
    zdamin=geq_2pi/(10*nbout-1)

    allocate(r1(inb1+1),z1(inb1+1))
    allocate(r2(inb2+1),z2(inb2+1))

    r1(1:inb1)=rb1; r1(inb1+1)=0.0_geq_r8
    z1(1:inb1)=zb1; z1(inb1+1)=0.0_geq_r8

    r2(1:inb2)=rb2; r2(inb2+1)=0.0_geq_r8
    z2(1:inb2)=zb2; z2(inb2+1)=0.0_geq_r8

    n1=inb1
    n2=inb2

    !  make sure contours are closed with duplicate points removed...

    call geq_ckrzc(r1,z1,n1,inb1+1)
    call geq_ckrzc(r2,z2,n2,inb2+1)

    rmaxwk=maxval(r1(1:n1))
    rminwk=minval(r1(1:n1))

    rcen1=0.55_geq_r8*rmaxwk + 0.45_geq_r8*rminwk  ! shift out 10%
    zcen1=0.5_geq_r8*(maxval(z1(1:n1))+minval(z1(1:n1)))

    rmaxwk=maxval(r2(1:n2))
    rminwk=minval(r2(1:n2))

    rcen2=0.55_geq_r8*rmaxwk + 0.45_geq_r8*rminwk  ! shift out 10%
    zcen2=0.5_geq_r8*(maxval(z2(1:n2))+minval(z2(1:n2)))

    rcen = zf1*rcen1 + zf2*rcen2
    zcen = zf1*zcen1 + zf2*zcen2

    allocate(a1(n1),d1(n1),a2(n2),d2(n2))

    call geq_aseq(rcen1,zcen1,r1,z1,n1,zdamin,a1,d1,n1a)
    call geq_aseq(rcen2,zcen2,r2,z2,n2,zdamin,a2,d2,n2a)

    do i=1,nbout
       ztha=geq_2pi/2 - (i-1)*geq_2pi/(nbout-1)  ! pi ... 0 ... -pi
       call geq_afind(a1,d1,n1a,ztha,dfind1)
       call geq_afind(a2,d2,n2a,ztha,dfind2)
       dfind = zf1*dfind1 + zf2*dfind2

       rbout(i)=rcen+cos(ztha)*dfind
       zbout(i)=zcen+sin(ztha)*dfind
    enddo

  end subroutine geq_mergbb

  subroutine geq_aseq(rcen,zcen,r,z,n,damin,a,d,na)

    ! convert closed sequence to succession of angles and distances
    ! relative to indicated center; force minimum angle change per step

    real(geq_r8), intent(in) :: rcen,zcen
    integer, intent(in) :: n
    real(geq_r8), intent(in) :: r(n),z(n)
    real(geq_r8), intent(in) :: damin
    real(geq_r8) :: a(n),d(n)   ! angles and distances
    integer, intent(out) :: na  ! number spanning +/- 2pi range

    !----------------------
    integer :: isign,i,ishift
    real(geq_r8) :: alast,alim,aa,aprev
    !----------------------

    isign=0
    do i=1,n-1
       a(i)=atan2(z(i)-zcen,r(i)-rcen)
       d(i)=sqrt((z(i)-zcen)**2+(r(i)-rcen)**2)
       if(i.gt.1) then
          if(a(i).gt.a(i-1)) then
             isign=isign+1
          else
             isign=isign-1
          endif
       endif
    enddo

    if(isign.lt.0) then
       isign=-1
       a(1:n-1)=a(1:n-1)*isign
    else
       isign=1
    endif
    
    alast=a(1)+geq_2pi
    alim=alast-damin
    a(n)=alast
    d(n)=d(1)

    ishift=0
    na=1
    i=1
    aprev=a(1)
    do
       i=i+1
       if(i.eq.n) then
          na=na+1
          a(na)=alast
          d(na)=d(n)
          exit
       endif
       aa=a(i)+ishift*geq_2pi
       if(aa.lt.aprev-geq_2pi/2) then
          ishift=1
          aa=aa+geq_2pi
       endif
       if((aa.ge.aprev+damin).and.(aa.le.alim)) then
          na=na+1
          a(na)=aa
          d(na)=d(i)
       endif
    enddo

    if(na.lt.n) then
       a(na+1:n)=0
       d(na+1:n)=0
    endif

    if(isign.eq.-1) a=a*isign

  end subroutine geq_aseq

  subroutine geq_afind(a,d,n,ang,dist)
    ! interpolate distance-- linear interpolation of parabolic fits
    !  as a function of angle

    integer, intent(in) :: n
    real(geq_r8), intent(in) :: a(n),d(n)  ! dist vs. angle
    real(geq_r8), intent(in) :: ang        ! angle in
    real(geq_r8), intent(out) :: dist      ! distance out

    ! a is either monotonic increasing or monotonic decreasing
    !----------------------------
    real(geq_r8) :: zang,zfac1,zfac2,zd1,zd2
    integer :: isign,iest
    !----------------------------
    ! force angle into range

    zang=ang

    if(zang.lt.min(a(1),a(n))) then
       do
          zang = zang + geq_2pi
          if(zang.ge.min(a(1),a(n))) exit
       enddo
    else if(zang.gt.max(a(1),a(n))) then
       do
          zang = zang - geq_2pi
          if(zang.le.max(a(1),a(n))) exit
       enddo
    endif

    if(a(n).gt.a(1)) then
       isign=1
    else
       isign=-1
    endif

    iest = 1 + (zang-a(1))*(n-1)/(a(n)-a(1))
    iest = max(1,min((n-1),iest))

    if(isign*zang.lt.isign*a(iest)) then
       do
          if(iest.eq.1) exit
          iest=iest-1
          if(isign*zang.ge.isign*a(iest)) exit
       enddo
    else if(isign*zang.gt.isign*a(iest+1)) then
       do
          if((iest+1).eq.n) exit
          iest=iest+1
          if(isign*zang.le.isign*a(iest+1)) exit
       enddo
    endif

    zfac2=(zang-a(iest))/(a(iest+1)-a(iest))
    zfac1=1.0_geq_r8-zfac2

    if(iest.eq.1) then
       call dfit(a(n-1)-isign*geq_2pi,d(n-1),a(1),d(1),a(2),d(2), &
            zang,zd1)
    else
       call dfit(a(iest-1),d(iest-1),a(iest),d(iest),a(iest+1),d(iest+1), &
            zang,zd1)
    endif

    if(iest.eq.(n-1)) then
       call dfit(a(n-1),d(n-1),a(n),d(n),a(2)+isign*geq_2pi,d(2), &
            zang,zd2)
    else
       call dfit(a(iest),d(iest),a(iest+1),d(iest+1),a(iest+2),d(iest+2), &
            zang,zd2)
    endif

    dist = zfac1*zd1 + zfac2*zd2

    contains
      subroutine dfit(am,dm,a,d,ap,dp,ang,dist)
        ! parabolic fit; ang btw am and ap; am,a,ap either strict increasing
        ! or strict decreasing.
        !   d(am)=dm; d(a)=d; d(ap)=dp; want d(ang)

        real(geq_r8), intent(in) :: am,dm,a,d,ap,dp,ang
        real(geq_r8), intent(out) :: dist

        real(geq_r8) :: xm,xp,fm,fp,c1,c2,x,den,f
        real(geq_r8) :: dlin,dmin,dmax

        !  f = c1*x**2 + c2*x;  ang = x+a;  dist = f+d

        xm = am-a
        xp = ap-a

        fm = dm-d
        fp = dp-d

        x = ang-a

        den = xp*xm*(xp-xm)
        c1 = (fp*xm-fm*xp)/den
        c2 = (fm*xp*xp-fp*xm*xm)/den

        f = x*(c1*x+c2)

        dist = f + d

        ! mod DMC: limits need to be imposed, for some extreme cases...
        ! these are based on the linear interpolation

        if(x*xm.ge.0.0_geq_r8) then
           f=x/xm
           dlin = f*dm + (1.0_geq_r8-f)*d
        else
           f=x/xp
           dlin = f*dp + (1.0_geq_r8-f)*d
        endif

        dmin = dlin*0.90_geq_r8
        dmax = dlin*1.15_geq_r8

        if((dist.lt.dmin).OR.(dist.gt.dmax)) then
           ! limits...
           !!! write(6,*) 'dm,d,dp = ',dm,d,dp
           !!! write(6,*) ' --> got: ',dist
           dist=max(dmin,min(dmax,dist))
           !!! write(6,*) ' --> use: ',dist
           !!! write(6,*) ' x_ratio: ',min(abs(xm),abs(xp))/max(abs(xm),abs(xp))
           !!! write(6,*) ' '
        endif

      end subroutine dfit

  end subroutine geq_afind

!------------------------------------------------------
  subroutine geq_gtime_num(geq,nt)

    ! fetch size of time vector

    type(geqdsk), intent(in) :: geq ! geqdsk object w/possible MDS+ time series
    integer, intent(out) :: nt      ! #of known times

    if(.not.allocated(geq%gtime_s)) then
       nt=0
    else
       nt=geq%t_num_valid
    endif

  end subroutine geq_gtime_num

  subroutine geq_gtime(geq,ztime,idim,nt)

    ! fetch time vector

    type(geqdsk), intent(in) :: geq ! geqdsk object w/possible MDS+ time series
    integer, intent(in) :: idim     ! dimension of time array
    real(geq_r8), intent(out) :: ztime(idim)      ! time array
    integer, intent(out) :: nt      ! #of known times

    ztime=0

    if(.not.allocated(geq%gtime_s)) then
       nt=0
    else
       nt=geq%t_num_valid
       if(idim.ge.nt) then
          ztime(1:nt)=geq%gtime_s(1:nt)
       endif
    endif

  end subroutine geq_gtime

  subroutine geq_curtime(geq,ztime1)

    ! fetch the current time (if MDS+ EFIT tree is in use)
    !  if data was from a G file, no time is available; ztime1=0.00 will be
    !  returned.

    type(geqdsk), intent(in) :: geq ! geqdsk object w/possible MDS+ time series
    real(geq_r8), intent(out) :: ztime1  ! time returned

    if(.not.allocated(geq%gtime_s)) then
       ztime1=0.00
    else
       ztime1=geq%curtime
    endif

  end subroutine geq_curtime

!------------------------------------------------------
!  parse MDS+ id header; look for "MDS+" leader, option name,
!   and keywords REDUCE or NOREDUCE
!
!   examples of valid header string:
!
!      MDS+
!      MDS+,REDUCE
!      MDS+,NOREDUCE
!      MDS+,NSTX
!      MDS+,NSTX,NOREDUCE
!
  subroutine geq_prepars(zbuf,opt,ireduce,ier)
    implicit NONE
    character*(*), intent(in) :: zbuf    ! input headerstring
    character*(*), intent(out) :: opt    ! option name
    logical, intent(out) :: ireduce      ! .TRUE. = "REDUCE"
    integer, intent(out) :: ier          ! status code, 0=OK
 
! local...
 
    integer idlima,idlim2,ilen,str_length
    integer iw1,iw2
    character*10 ztest,zword1,zword2
 
! executable code...
 
    ier=0
    ireduce=.TRUE.         ! default:  use REDUCE option
    opt=' '                ! default:  no option specified
 
    ztest=zbuf(1:4)
    call uupper(ztest)
    if(ztest.ne.'MDS+') then
       ier=3
       return              ! must start with "MDS+"
    endif
 
    ilen=str_length(zbuf)
    if(ilen.eq.4) return   ! bare "MDS+"
 
    idlima=index(zbuf,'/')
    if((idlima.eq.0).or.(idlima.eq.ilen)) then
       ier=3               ! slash can't be the last character;
       return              ! need slashes if more than just "MDS+"
    endif
 
    call trmds_pshrink(zbuf,1,idlima-1,iw1,iw2,ier)
    if(ier.ne.0) return
    ztest=zbuf(iw1:iw2)
    call uupper(ztest)
    if(ztest.ne.'MDS+') then
       ier=3
       return              ! apparent extraneous characters after "MDS+"
    endif
!
!  OK there are one or two words after the "MDS+"
!
    idlim2=index(zbuf(idlima+1:ilen),'/')
    if(idlim2.eq.0) then
!  one word only
       call trmds_pshrink(zbuf,idlima+1,ilen,iw1,iw2,ier)
       if(ier.ne.0) return
       zword1=zbuf(iw1:iw2)
       call uupper(zword1)
       zword2=' '
    else
!  two words
       idlim2=idlima+idlim2
       call trmds_pshrink(zbuf,idlima+1,idlim2-1,iw1,iw2,ier)
       if(ier.ne.0) return
       zword1=zbuf(iw1:iw2)
       call uupper(zword1)
       call trmds_pshrink(zbuf,idlim2+1,ilen,iw1,iw2,ier)
       if(ier.ne.0) return
       zword2=zbuf(iw1:iw2)
       call uupper(zword2)
    endif
!
!  check for REDUCE/NOREDUCE keywords
!
    if(zword1.eq.'REDUCE') then
       ireduce=.TRUE.
       zword1=' '
    endif
    if(zword1.eq.'NOREDUCE') then
       ireduce=.FALSE.
       zword1=' '
    endif
    if(zword2.eq.'REDUCE') then
       ireduce=.TRUE.
       zword2=' '
    endif
    if(zword2.eq.'NOREDUCE') then
       ireduce=.FALSE.
       zword2=' '
    endif
!
!  at most one non-blank can remain
!
    if((zword1.eq.' ').and.(zword2.eq.' ')) return  ! done, OK, opt=' '
!
    if((zword1.ne.' ').and.(zword2.ne.' ')) then
       ier=3
       return             ! error:  two non-blank non-keywords
    endif
!
!  OK set the option string
!
    if(zword1.eq.' ') then
       opt=zword2
    else
       opt=zword1
    endif
!
    return
  end subroutine geq_prepars
 
!-----------------------------------------------
!  free memory associated with G-EQdisk object
!  if the object was never initialized: deallocate errors are ignored!
!
  subroutine geq_free(geq, ier)
    implicit none
    type(geqdsk), intent(inout) :: geq
    integer, intent(out) :: ier
 
    ier = 0
 
! basic G-EQdisk items
 
    call geq_free_slice(geq)
 
! items for time dependent G-EQdisk data from MDSplus
 
    call geq_free_tvar(geq)
 
! clear identity information
 
    geq%access_method = 0
    geq%ident = ' '
    geq%mdstime = 0
    geq%mdsopt = ' '
    geq%mdslice = .false.
    geq%mdstok = ' '
 
  end subroutine geq_free
 
!-----------------------------------------------
! deallocate G-EQdisk time slice data
 
  subroutine geq_free_slice(geq)
    implicit none
    type(geqdsk), intent(inout) :: geq
 
    integer istat
 
    deallocate(geq%fpol_tm,stat=istat)
    deallocate(geq%pres_Nt__M2,stat=istat)
    deallocate(geq%ffprim_T2M2Rad__Wb,stat=istat)
    deallocate(geq%pprime_NtRad__M2Wb,stat=istat)
    deallocate(geq%qpsi,stat=istat)
 
    deallocate(geq%psirz_Wb__Rad,stat=istat)
 
    deallocate(geq%rbbbs_m,stat=istat)
    deallocate(geq%zbbbs_m,stat=istat)
 
  end subroutine geq_free_slice
 
!-----------------------------------------------
! deallocate (MDS+) time variation data
 
  subroutine geq_free_tvar(geq)
    implicit none
    type(geqdsk), intent(inout) :: geq
 
    integer istat
 
    deallocate(geq%gtime_s,stat=istat)
 
    deallocate(geq%t_index,stat=istat)
    deallocate(geq%t_nw,stat=istat)
    deallocate(geq%t_nh,stat=istat)
 
    deallocate(geq%t_rdim,stat=istat)
    deallocate(geq%t_zdim,stat=istat)
    deallocate(geq%t_rcentr,stat=istat)
    deallocate(geq%t_rleft,stat=istat)
    deallocate(geq%t_zmid,stat=istat)
    deallocate(geq%t_rmaxis,stat=istat)
    deallocate(geq%t_zmaxis,stat=istat)
    deallocate(geq%t_bcentr,stat=istat)
    deallocate(geq%t_current,stat=istat)
    deallocate(geq%t_simag,stat=istat)
    deallocate(geq%t_sibry,stat=istat)
 
    deallocate(geq%t_fpol,stat=istat)
    deallocate(geq%t_pres,stat=istat)
    deallocate(geq%t_ffprim,stat=istat)
    deallocate(geq%t_pprime,stat=istat)
    deallocate(geq%t_qpsi,stat=istat)
 
    deallocate(geq%t_psirz,stat=istat)
 
    deallocate(geq%t_nbbbs,stat=istat)
 
    deallocate(geq%t_rbbbs,stat=istat)
    deallocate(geq%t_zbbbs,stat=istat)
 
    deallocate(geq%rlim_m,stat=istat)
    deallocate(geq%zlim_m,stat=istat)
 
  end subroutine geq_free_tvar
 
!.............................................................................
 
  subroutine geq_error(ier)
    implicit none
    integer, intent(in) :: ier
 
    integer lunout
 
!.....................
 
    call geqdsk_lunmsg_get(lunout)
 
    select case (ier)
 
       case (1)
          write(lunout,*) &
               & 'geqdsk::geqdsk_init:: ERROR cannot open G-EQDSK file!'
       case (2)
          write(lunout,*) &
               & 'geqdsk::geqdsk_init:: ERROR occurred while allocating!'
       case (3)
          write(lunout,*) &
               & 'geqdsk::geqdsk_init:: MDS+ access string parse error!'
       case (4)
          write(lunout,*) &
               & 'geqdsk::geqdsk_init:: MDS+ slice time request not readable!'
       case (5)
          write(lunout,*) &
               & 'geqdsk::geqdsk_init:: MDS+ error (unspecified)!'
       case (6)
          write(lunout,*) &
               & 'geqdsk::geqdsk_init:: MDS+ i/o error!'
       case (7)
          write(lunout,*) &
               & 'geqdsk::geqdsk_init:: MDS+ data consistency error!'
 
       case default
 
    end select
  end subroutine geq_error
 
!-----------------------------------------------
  subroutine geq_echo(str1,str2)
    implicit none
    character*(*) str1,str2   ! message components
 
    integer lunout,ilen,str_length
 
    call geqdsk_lunmsg_get(lunout)
 
    ilen=max(1,str_length(str1))
    write(lunout,'(a)') str1(1:ilen)
    ilen=max(1,str_length(str2))
    write(lunout,'(2x,a)') str2(1:ilen)
 
  end subroutine geq_echo
 
  subroutine geq_nanscan(zexpr,inum,zr4buf,ilnanflag)
    character*(*), intent(in) :: zexpr  ! MDSplus expression (for message)
    integer, intent(in) :: inum         ! array dimension
    real, intent(inout) :: zr4buf(inum)       ! array to be scanned for NANs
    logical, intent(inout) :: ilnanflag(inum) ! FLAG any NaN occurrences

    ! detect NaNs by writing value into a character string; if character string
    ! using "*" format.  If this contains no decimal point assume a NaN

    ! a warning message is printed, the element flag is set, and the value
    ! is cleared to prevent FPEs occurring later...

    ! NOTE that the ilnanflag array is NOT cleared to .FALSE. This is the
    ! caller's responsibility; this allows multiple geq_nanscan calls on
    ! multiple equivalently shaped data items, s.t. the flag is a logical OR
    ! of NaNs occurring in any of the data items.

    character*30 :: nantest
    integer :: ii,lunout

    do ii=1,inum
       nantest=' '
       write(nantest,*) zr4buf(ii)
       if(index(nantest,'.').eq.0) then
          ilnanflag(ii)=.TRUE.
          zr4buf(ii)=0.0
          call geqdsk_lunmsg_get(lunout)
          write(lunout,*) ' *** at index ',ii,': ',trim(nantest), &
               ' detected in ',trim(zexpr)
       endif
    enddo

  end subroutine geq_nanscan

  subroutine geq_bchk(geq,zlbl,inum,rb,zb,ier)

    ! scan the bdy data -- report error on any points "outside the box".

    use geqdsk_aux
    implicit NONE

    type(geqdsk), intent(in) :: geq    ! G-eqdsk data object
    character*(*), intent(in) :: zlbl  ! label for contour
    integer, intent(in) :: inum        ! number of (R,Z) points
    real*8, intent(inout) :: rb(inum),zb(inum)  ! the (R,Z) point sequence

    integer, intent(out) :: ier        ! return code, 0=OK

    real*8, parameter :: threshold = 2.d-5  ! fix up out of bound points which are less than this threshold

    !--------------------------
    ! local:
    real*8 :: zrmin,zrmax,zzmin,zzmax
    integer :: ib,iwarnR,iwarnZ,iwasR,iwasZ
    !--------------------------

    zrmin = geq%Rleft_m
    zrmax = zrmin + geq%Rdim_m

    zzmin = geq%Zmid_m - geq%Zdim_m * 0.5d0
    zzmax = zzmin + geq%Zdim_m

    ier = 0
    iwarnR = 0
    iwarnZ = 0

    do ib = 1,inum
       if((rb(ib).lt.zrmin).or.(rb(ib).gt.zrmax)) then
          iwasR = iwarnR
          if(iwarnR.eq.0) then
             iwarnR=1
             write(lunmsg,*) ' ?geq_bchk: '//trim(zlbl)//' out of R bounds:'
             write(lunmsg,*) '  Rmin, Rmax = ',zRmin,zRmax
          endif
          write(lunmsg,*) '  Rb(',ib,')= ',Rb(ib)
          if (iwasR<=0) then
             if (max(zrmin-rb(ib),rb(ib)-zrmax)<threshold*(zrmax-zrmin)) then
                write(lunmsg,*) '  ---> !fixup applied'
                rb(ib) = min(zrmax,max(rb(ib),zrmin))
                iwarnR = -1
             else
                write(lunmsg,*) '  ---> ?fixup not possible'
             end if
          end if
       endif
       if((zb(ib).lt.zzmin).or.(zb(ib).gt.zzmax)) then
          iwasZ = iwarnZ
          if(iwarnZ.eq.0) then
             iwarnZ=1
             write(lunmsg,*) ' ?geq_bchk: '//trim(zlbl)//' out of Z bounds:'
             write(lunmsg,*) '  Zmin, Zmax = ',zZmin,zZmax
          endif
          write(lunmsg,*) '  Zb(',ib,')= ',Zb(ib)
          if (iwasZ<=0) then
             if (max(zzmin-zb(ib),zb(ib)-zzmax)<threshold*(zzmax-zzmin)) then
                write(lunmsg,*) '  ---> !fixup applied'
                zb(ib) = min(zzmax,max(zb(ib),zzmin))
                iwarnZ = -1
             else
                write(lunmsg,*) '  ---> ?fixup not possible'
             end if
          end if
       endif
    enddo

    if (iwarnR>0 .or. iwarnZ>0) ier=1
  end subroutine geq_bchk

end module geqdsk_mod
