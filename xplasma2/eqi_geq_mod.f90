module eqi_geq_mod
 
  use ezspline_obj
  use ezspline
  use geqdsk_mod
 
  !  DMC 2000 -- 2008 ... ad hoc data for contouring EFIT-style Psi(R,Z) data

  implicit NONE
  save  !  caution -- not saved/restored by eq_save/eq_restore
 
  type(ezspline2_r8) :: pspl   ! psi(R,Z) spline
 
  type(geqdsk) :: geq          ! GEQDISK data
 
  real*8 psep                  ! EFIT's Psi @ bdy

  real*8 r0,z0                 ! mag. axis
 
  real*8 :: gs_errmax          ! maximum GS error found in EFIT equilibrium
                               ! (ignoring pts near bdy)
 
  integer :: jwarn = 0         ! warning flag

  integer :: isign_xy          ! trace direction for contouring
  integer :: errflag           ! set if 

  real*8 :: lnorm,lmax         ! normalizing distance for contouring
                               ! or for curvature definition

  integer :: isw               ! toggle btw 2 th_adapt grids & assoc. data
  real*8, dimension(:,:), allocatable :: thadapt
                               ! adaptive th grid -- will concentrate points
                               ! near approach to bdy separatrices in Psi(R,Z)

  integer :: nthad             ! size of data in thadapt(:,isw)
  integer :: nthad0            ! starting size of data in thadapt(:,isw)

  real*8, dimension(:,:), allocatable :: radapt,zadapt,curvad,rzwk
  !  radapt(...) R pts @ thadapt on a chosen surface
  !  zadapt(...) Z pts @ thadapt on a chosen surface
  !  curvad(...) normalized contour curvature on the chosen surface
  !  rzwk(...) work space for eqi_find_curvature routine (see also curvec(...))

  integer, dimension(:), allocatable :: iadapt

  !  auxilliary utility (R,Z) countour
  integer :: nseek
  real*8, dimension(:), allocatable :: rseek,zseek,curvec

  real*8 :: cbdy_chk

  integer :: istep_out

contains
  subroutine nx_thadapt(inref,inmin)

    integer, intent(in) :: inref
    !  inref -- reference gridsize -- faster adaptation if at or below...
    integer, intent(in) :: inmin
    !  inmin -- minimum size of contigous region with curvature <= cbdy_chk
    !           lesser size => grid adaptation

    !----------------------------
    integer :: iswp,imul,ii,jj,kk,ict,indx1,indxn,indxi,indxo,iaprev
    real*8 :: mulfac
    integer, save :: icount = 0
    integer, parameter :: ilun_dbg = 6
    !----------------------------

    icount = icount + 1

    iswp=isw
    isw=3-isw

    if(nthad.le.inref) then
       imul=10
    else
       imul=2
    endif
    mulfac = 1.0d0
    mulfac = mulfac/imul

    iaprev=1
    do ii=1,nthad
       iadapt(ii)=0
       if(curvad(ii,iswp).le.cbdy_chk) iadapt(ii)=1
       if((iadapt(ii).eq.1).and.(iaprev.eq.0)) then
          iadapt(ii-1)=1 ! extend 1 back
       endif
       if(ii.eq.1) iaprev=0
       if((iadapt(ii).eq.0).and.(iaprev.eq.1)) then
          iaprev=0
          iadapt(ii)=1   ! extend forward
       else
          iaprev=iadapt(ii)
       endif
    enddo

    !  see if adapt region straddles end pt
    indx1=0
    indxn=nthad

    if(iadapt(1)*iadapt(nthad).eq.1) then
       ict=-1   ! do not double count common end point
       ii = 1
       do while ( iadapt(ii) == 1 )
          ict = ict + 1
          ii = ii + 1
       enddo
       indx1 = ii-1

       ii = nthad
       do while ( iadapt(ii) == 1 )
          ict = ict + 1
          ii = ii - 1
       enddo
       indxn = ii + 1
    endif
    if(indx1.gt.0) then
       if(ict.ge.inmin) then
          indx1=0
          indxn=nthad
       endif
    endif

    indxi=0
    indxo=0
    ! adapt 1st piece contigous with index 1 if necessary
    do ii=1,indx1
       indxi=indxi+1
       indxo=indxo+1
       thadapt(indxo,isw)=thadapt(indxi,iswp)
       do jj=1,imul-1
          indxo=indxo+1
          thadapt(indxo,isw)=((imul-jj)*thadapt(indxi,iswp) + &
               jj*thadapt(indxi+1,iswp))*mulfac
       enddo
    enddo

    ii=indx1
    ict = 0
    do while ( ii < indxn-1 )
       ! look for more contiguous sections with low curvature radius
       ii = ii + 1
       if(iadapt(ii).eq.0) then
          if(ict.eq.0) then
             ! no adaptation on this point...
             continue
          else if(ict.ge.inmin) then
             ! contiguous section large enough; no adaptation
             do jj=1,ict
                indxi=indxi+1
                indxo=indxo+1
                thadapt(indxo,isw)=thadapt(indxi,iswp)
             enddo
          else
             ! contiguous low curvature region has too few points -- expand
             do jj=1,ict
                indxi=indxi+1
                indxo=indxo+1
                thadapt(indxo,isw)=thadapt(indxi,iswp)
                do kk=1,imul-1
                   indxo=indxo+1
                   thadapt(indxo,isw)=((imul-kk)*thadapt(indxi,iswp) + &
                        kk*thadapt(indxi+1,iswp))*mulfac
                enddo
             enddo
          endif
          ict = 0  ! clear count
          ! include current point
          indxi=indxi+1
          indxo=indxo+1
          thadapt(indxo,isw)=thadapt(indxi,iswp)
       else
          ict = ict + 1  ! keep count
       endif
    enddo

    ! adapt last piece contiguous with index 
    do ii=indxn,nthad-1
       indxi=indxi+1
       indxo=indxo+1
       thadapt(indxo,isw)=thadapt(indxi,iswp)
       do jj=1,imul-1
          indxo=indxo+1
          thadapt(indxo,isw)=((imul-jj)*thadapt(indxi,iswp) + &
               jj*thadapt(indxi+1,iswp))*mulfac
       enddo
    enddo

    ! final end point
    indxi=indxi+1
    if(indxi.ne.nthad) then
       write(ilun_dbg,*) &
            '    nx_thadapt: Error occurred on call icount = ',icount
       call errmsg_exit(' ?? eqi_geq_mod(nx_thadapt) algorithm failure.')
    endif
    indxo=indxo+1

    thadapt(indxo,isw)=thadapt(indxi,iswp)
    nthad=indxo

  end subroutine nx_thadapt

  subroutine quad_minmax(delx12,delx23,f1,f2,f3,fminmax)

    ! given a sequence of 3 consecutive function values {f1,f2,f3} evaluated
    ! at grid points x1 < x2 < x3 (with delx12 = x2-x1; delx23 = x3-x2), and
    ! with f2 a pointwise local minimum or maximum, find the local minimum
    ! or maximum based on a parabolic fit through the 3 points.

    !    f(-delx1) = f2 + aa*delx1**2 - bb*delx1 = f1
    !    f(delx2)  = f2 + aa*delx2**2 + bb*delx2 = f3

    ! f2 **must** be either .le.min(f1,f3) or .ge.max(f1,f3) -- not checked!!

    real*8, intent(in) :: f1,f2,f3,delx12,delx23
    real*8, intent(out) :: fminmax

    real*8 :: xm,xp,fm,fp,aa,bb,denom,xsing
    real*8, parameter :: zero = 0.0d0

    xm = -delx12
    xp =  delx23
    fm = f1-f2
    fp = f3-f2

    denom=xp*xm*(xp-xm)
    aa=(fp*xm-fm*xp)/denom
    bb=(fm*xp*xp-fp*xm*xm)/denom

    if(aa.eq.zero) then
       fminmax=f2   ! assume flat, colinear points f1=f2=f3
    else
       xsing = -bb/(2*aa)
       fminmax = f2 + xsing*(bb + xsing*aa)
    endif

  end subroutine quad_minmax

end module eqi_geq_mod
