subroutine xplasma_integ_work(integ_num,nth,thvec,wth, &
     nths,ths,nrhos,rhos,rhomin, &
     iavail,wk,dvdrho,bmaxa,result,ier)

  !  compute surface oriented numeric integration, e.g. flux surface 
  !  integrals and averages to yield results of form f(rho) or f(theta,rho)

  !  because certain terms (R, dR/dtheta, etc.) come up repeatedly, a
  !  caching mechanism is supported.

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  !-----------------------------------------------------------------------
  integer, intent(in) :: integ_num   ! id of integral/average to be done

  integer, intent(in) :: nth                ! #of theta evaluations
  real*8, intent(in) :: thvec(nth)          ! theta points
  real*8, intent(in) :: wth(nth)            ! associated d(theta)*weights

  integer, intent(in) :: nths               ! #of theta zone bdys
  real*8, intent(in) :: ths(nths)           ! theta zone boundaries
  integer, intent(in) :: nrhos              ! #of rho surfaces
  real*8, intent(in) :: rhos(nrhos)         ! rho surfaces
  real*8, intent(in) :: rhomin              ! minimum safe rho for evaluation

  integer, parameter :: ifuns=10            ! #of functions in cache

  logical, intent(inout) :: iavail(ifuns)      ! cache status flags
  real*8, intent(inout) :: wk(nth,ifuns,nrhos) ! cache

  real*8, intent(in) :: dvdrho(nths-1,nrhos)   ! dV/drho at surfaces (for avgs)
  real*8, intent(inout) :: bmaxa(nrhos)        ! max(|B|) on surfaces
  !  if this is zero, and it is needed here, it will be evaluated...

  real*8, intent(out) :: result(nths-1,nrhos)  ! result of integration or avg
  integer, intent(out) :: ier         ! status code, 0=OK

  !  avg is short for "average".
  !-----------------------------------------------------------------------
  !  local...

  type (xpeval) :: th_intrp, rho_intrp

  integer :: iertmp
  logical :: lav_dvdrho,iflag10,local_eval
  integer :: id_R,id_Z,id_Bmod,id_BR,id_BZ,id_rho,id_th
  character*100 zmsg

  ! cache line codes:
  !   1-- R  2-- dR/dtheta  3-- dR/drho
  !   4-- Z  5-- dZ/dtheta  6-- dZ/drho
  !   7-- Bmod   8-- BR     9-- BZ
  !       10-- R*|det(J)| = 
  !              2*pi*R*abs((dR/drho)*(dZ/dtheta)-(dZ/drho)*(dR/dtheta))
  integer :: icad(ifuns),inca,ica,ii,ith,ithr,ieval,irhos

  integer :: inth0,indexth(nth)
  real*8 :: thdum(nth/10),rhodum(nth/10),fintdum(nth/10)
  real*8, dimension(:), allocatable :: fint, local_wk
  real*8, parameter :: c2pi = 6.2831853071795862D+00

  integer, save :: integ_count = 0  ! for debugging...
  integer, parameter :: debug_count1 = -1
  integer, parameter :: debug_count2 = -1
  integer, parameter :: ilun_dbg=6
  !-----------------------------------------------------------------------

  ier=0

  call xplasma_common_ids(sp,ier, id_R=id_R, id_Z=id_Z, &
       id_Bmod=id_Bmod, id_BR=id_BR, id_BZ=id_BZ)

  if(ier.ne.0) return

  call xplasma_prof_info(sp, id_R, ier,  gridId1=id_th, gridId2=id_rho)

  if(ier.ne.0) return

  call xplasma_integ_mkindex(nths,ths,nth,thvec,indexth)

  ! for averages at the axial singularity...
  inth0 = nth/10
  do ith=1,inth0
     thdum(ith)=(ith-1)*c2pi/inth0
     rhodum(ith)=rhomin
  enddo

  if(id_Bmod.eq.0) call xplasma_find_item(sp,"BMOD",id_Bmod,iertmp,nf_noerr=.TRUE.)
  if(id_BR.eq.0) call xplasma_find_item(sp,"BR",id_BR,iertmp,nf_noerr=.TRUE.)
  if(id_BZ.eq.0) call xplasma_find_item(sp,"BZ",id_BZ,iertmp,nf_noerr=.TRUE.)

  iflag10=.FALSE.
  local_eval=.FALSE.

  ! determine whether an integral or an average; determine what data
  ! items (R,Z, derivatives of R,Z, mod(B), B components...) are needed.

  if(integ_num.eq.eqi_intg_vol) then
     lav_dvdrho=.FALSE.
     inca=3
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_z

  else if(integ_num.eq.eqi_intg_dvdrho) then
     lav_dvdrho=.FALSE.
     inca=5
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_drdrho
     icad(4)=eqi_intg_dzdth; icad(5)=eqi_intg_dzdrho
     iflag10=.TRUE.

  else if((integ_num.eq.eqi_intg_area).or.(integ_num.eq.eqi_intg_darea)) then
     lav_dvdrho=.FALSE.
     inca=2
     icad(1)=eqi_intg_drdth; icad(2)=eqi_intg_z

  else if(integ_num.eq.eqi_intg_r2i) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_r

  else if(integ_num.eq.eqi_intg_itor) then
     lav_dvdrho=.FALSE.
     inca=4
     icad(1)=eqi_intg_drdth; icad(2)=eqi_intg_dzdth
     icad(3)=eqi_intg_br; icad(4)=eqi_intg_bz

  else if(integ_num.eq.eqi_intg_lpol) then
     lav_dvdrho=.FALSE.
     inca=2
     icad(1)=eqi_intg_drdth; icad(2)=eqi_intg_dzdth

  else if(integ_num.eq.eqi_intg_surf) then
     lav_dvdrho=.FALSE.
     inca=3
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_dzdth

  else if(integ_num.eq.eqi_intg_booz_qovg) then
     lav_dvdrho=.FALSE.
     inca=5
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_dzdth
     icad(4)=eqi_intg_br; icad(5)=eqi_intg_bz

  else if(integ_num.eq.eqi_intg_g2r2i) then
     lav_dvdrho=.TRUE.
     inca=5
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_drdrho
     icad(4)=eqi_intg_dzdth; icad(5)=eqi_intg_dzdrho

  else if(integ_num.eq.eqi_intg_ri) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_r

  else if(integ_num.eq.eqi_intg_r3i) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_r

  else if(integ_num.eq.eqi_intg_r2) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_r

  else if(integ_num.eq.eqi_intg_r1) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_r

  else if(integ_num.eq.eqi_intg_br2) then
     lav_dvdrho=.TRUE.
     inca=2
     icad(1)=eqi_intg_r
     icad(2)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_bi) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_b2i) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_b2) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_b1) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_bz2) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_bz

  else if(integ_num.eq.eqi_intg_g2b2i) then
     lav_dvdrho=.TRUE.
     inca=5
     icad(1)=eqi_intg_drdth; icad(2)=eqi_intg_drdrho
     icad(3)=eqi_intg_dzdth; icad(4)=eqi_intg_dzdrho
     icad(5)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_g1) then
     lav_dvdrho=.TRUE.
     inca=4
     icad(1)=eqi_intg_drdth; icad(2)=eqi_intg_drdrho
     icad(3)=eqi_intg_dzdth; icad(4)=eqi_intg_dzdrho

  else if(integ_num.eq.eqi_intg_g2) then
     lav_dvdrho=.TRUE.
     inca=4
     icad(1)=eqi_intg_drdth; icad(2)=eqi_intg_drdrho
     icad(3)=eqi_intg_dzdth; icad(4)=eqi_intg_dzdrho

  else if(integ_num.eq.eqi_intg_g2r3i) then
     lav_dvdrho=.TRUE.
     inca=5
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_drdrho
     icad(4)=eqi_intg_dzdth; icad(5)=eqi_intg_dzdrho

  else if(integ_num.eq.eqi_intg_g2r2) then
     lav_dvdrho=.TRUE.
     inca=5
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_drdrho
     icad(4)=eqi_intg_dzdth; icad(5)=eqi_intg_dzdrho

  else if(integ_num.eq.eqi_intg_giri) then
     lav_dvdrho=.TRUE.
     inca=5
     icad(1)=eqi_intg_r; icad(2)=eqi_intg_drdth; icad(3)=eqi_intg_drdrho
     icad(4)=eqi_intg_dzdth; icad(5)=eqi_intg_dzdrho

  else if(integ_num.eq.eqi_intg_fhnclass) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_fhnc_tsc) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_modb

  else if(integ_num.eq.eqi_intg_fhmmx) then
     lav_dvdrho=.TRUE.
     inca=1
     icad(1)=eqi_intg_modb

  else if((integ_num.gt.eqi_intg_bign).and. &
       (integ_num.le.2*eqi_intg_bign)) then

     ! user defined average

     lav_dvdrho=.TRUE.
     inca=0
     local_eval=.TRUE.

  else if(integ_num.gt.2*eqi_intg_bign) then

     ! user defined integral

     lav_dvdrho=.FALSE.
     inca=0
     local_eval=.TRUE.

  else
     
     ier=611
     zmsg=' '
     write(zmsg,*) '   Passed integration code unrecognized in xplasma_integ_work: ',integ_num
     call xplasma_errmsg_append(sp,zmsg)
     return
  endif

  if(lav_dvdrho.and.(.not.iavail(10))) then
     ier=611
     call xplasma_errmsg_append(sp, &
         ' logic error in xplasma_integ_work: dVdrho integrand not available.')
     return
  endif

  allocate(fint(nth),local_wk(nth))

  ieval=0
  do ii=1,inca
     ica=icad(ii)
     if(.not.iavail(ica)) then
        ieval=ieval+1
     endif
  enddo

  if(local_eval.or.(ieval.gt.0)) then
     !  do theta lookups now -- these do not vary with (rho)

     call xplasma_x_lookup(sp, id_th, thvec, th_intrp, ier)
     if(ier.ne.0) then
        call xpeval_free(th_intrp)
        return
     endif
     
  endif

  do irhos=1,nrhos

     integ_count = integ_count + 1

     if(integ_count.eq.debug_count1) then
        write(ilun_dbg,*) ' reached debug_count1!'
     endif

     if(integ_count.eq.debug_count2) then
        write(ilun_dbg,*) ' reached debug_count2!'
     endif

     if(iflag10.and.iavail(10)) then
        continue                     ! needed data already in hand...

     else

        if(local_eval.or.(ieval.gt.0)) then
           call xplasma_x_lookup(sp, id_rho, rhos(irhos), rho_intrp, ier)
           if(ier.ne.0) then
              call xpeval_free(rho_intrp)
              call xpeval_free(th_intrp)
              return
           endif
        endif

        if(ieval.gt.0) then
           do ii=1,inca
              ica=icad(ii)
              if(iavail(ica)) cycle

              if(rhos(irhos).le.rhomin) then
                 call xplasma_integ_wk1(ica,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
                      inth0,thdum,rhodum,wk(1,ica,irhos),ier)
              else
                 call xplasma_integ_wk2(ica,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
                      nth,th_intrp,rho_intrp,wk(1,ica,irhos),ier)
              endif
              if(ier.ne.0) exit
           enddo
        endif
        if(ier.ne.0) exit

        if(local_eval) then
           if((rhos(irhos).le.rhomin).and.lav_dvdrho) then
              call xplasma_integ_wk1(integ_num,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
                   inth0,thdum,rhodum,local_wk,ier)
           else
              call xplasma_integ_wk2(integ_num,id_R,id_Z,id_Bmod,id_BR,id_BZ, &
                   nth,th_intrp,rho_intrp,local_wk,ier)
           endif
        endif
        if(ier.ne.0) exit

        result(1:nths-1,irhos)=0

        if(rhos(irhos).le.rhomin) then
           if(lav_dvdrho) then
              wk(1:inth0,10,irhos)=1 ! even weighting vs. theta
              call xplasma_integ_fwkdv(integ_num,inth0,nth,wk(1,1,irhos), &
                   local_wk,bmaxa(irhos),fintdum,ier)
              result(1:nths-1,irhos)=sum(fintdum)/inth0 ! avg axial value of
           else

           !  result stays zero for integrals (not averages) @rho=rhomin
           !  exception: 

              if(local_eval) then
                 call xplasma_integ_fwk(integ_num,nth,wk(1,1,irhos),local_wk, &
                      fint,ier)

                 do ith=1,nth
                    ithr=indexth(ith)
                    result(ithr,irhos)=result(ithr,irhos)+fint(ith)*wth(ith)
                 enddo
              endif

           endif

        else
           if(lav_dvdrho) then
              call xplasma_integ_fwkdv(integ_num,nth,nth,wk(1,1,irhos), &
                   local_wk,bmaxa(irhos),fint,ier)
           else
              call xplasma_integ_fwk(integ_num,nth,wk(1,1,irhos),local_wk, &
                   fint,ier)

              if(iflag10) then
                 wk(1:nth,10,irhos) = fint
              endif

           endif

           do ith=1,nth
              ithr=indexth(ith)
              result(ithr,irhos)=result(ithr,irhos)+fint(ith)*wth(ith)
           enddo

           if(lav_dvdrho) then
              do ithr=1,nths-1
                 result(ithr,irhos)=result(ithr,irhos)/dvdrho(ithr,irhos)
              enddo
           endif
        endif
     endif
     if(ier.ne.0) exit

  enddo

  if((ieval.gt.0).and.(ier.eq.0)) then
     do ii=1,inca
        ica=icad(ii)
        iavail(ica)=.TRUE.
     enddo
  endif

  if(iflag10.and.(ier.eq.0)) then
     iavail(10)=.TRUE.
  endif

  deallocate(fint,local_wk)

  call xpeval_free(rho_intrp)
  call xpeval_free(th_intrp)

end subroutine xplasma_integ_work

subroutine xplasma_integ_gwork(inum,inthbdys,thbdys, &
     inrhos,rhos,infine,rhofine,wrhof,indexs,th_result,result,ier)

  ! complete Guass-integration of collection of 2d zones' volumes or
  ! cross sectional areas

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  !-----------------------------------------------------------------------

  integer, intent(in) :: inum      ! integration code number

  integer, intent(in) :: inthbdys  ! theta boundaries of zones
  real*8, intent(in) :: thbdys(inthbdys)

  integer, intent(in) :: inrhos    ! rho boundaries of zones
  real*8, intent(in) :: rhos(inrhos)

  integer, intent(in) :: infine    ! rho integration grid, weights, rhos index
  real*8, intent(in) :: rhofine(infine)
  real*8, intent(in) :: wrhof(infine)
  integer, intent(in) ::  indexs(infine)

  !  result of partial integration in theta direction
  real*8, intent(in) :: th_result(inthbdys-1,inrhos)

  !  completed result
  real*8, intent(out) :: result(inthbdys-1,inrhos-1)

  integer, intent(out) :: ier           ! completion code,  0=OK

  !---------------------
  real*8, dimension(:,:), allocatable :: zwk
  real*8, dimension(:), allocatable :: rr,zz,drdrho
  integer :: irho,ii,ith,id_R,id_Z,id_th,id_rho
  real*8 :: zincr

  type (xpeval) :: th_intrp, rho_intrp
  real*8, parameter :: c2pi = 6.2831853071795862D+00
  !---------------------
  
  ier=0

  call xplasma_common_ids(sp,ier, id_R=id_R, id_Z=id_Z)
  if(ier.ne.0) return

  call xplasma_prof_info(sp, id_R, ier,  gridId1=id_th, gridId2=id_rho)

  if((inum.ne.eqi_intg_darea).and.(inum.ne.eqi_intg_dvol)) then
     ier=9999
     call xplasma_errmsg_append(sp, ' xplasma_integ_gwork detected error:')
     call xplasma_errmsg_append(sp, ' dArea or dVol not being computed.')
     return
  endif

  allocate(zwk(inthbdys,inrhos-1),rr(infine),zz(infine),drdrho(infine))

  ! one rho lookup for all evaluations...

  call xplasma_x_lookup(sp, id_rho, rhofine, rho_intrp, ier)
  if(ier.ne.0) then
     call xpeval_free(rho_intrp)
     return
  endif

  do ith=1,inthbdys-1

     call xplasma_x_lookup(sp, id_th, thbdys(ith), th_intrp, ier)
     if(ier.ne.0) exit

     if(inum.ne.eqi_intg_darea) then
        call xplasma_eval_prof(sp,id_R, th_intrp, rho_intrp, rr, ier)
        if(ier.ne.0) exit
     endif

     call xplasma_eval_prof(sp,id_R, th_intrp, rho_intrp, drdrho, ier, &
          ideriv2=1)  ! dR/drho
     if(ier.ne.0) exit

     call xplasma_eval_prof(sp,id_Z, th_intrp, rho_intrp, zz, ier)
     if(ier.ne.0) exit

     zwk(ith,1:inrhos-1)=0
     do ii=1,infine
        if(inum.eq.eqi_intg_darea) then
           ! dArea contribution
           zincr = zz(ii)*drdrho(ii)*wrhof(ii)

        else if(inum.eq.eqi_intg_dvol) then
           ! dVol contribution
           zincr = c2pi*rr(ii)*zz(ii)*drdrho(ii)*wrhof(ii)
        endif

        irho=indexs(ii)
        zwk(ith,irho)=zwk(ith,irho)+zincr
     enddo

     if(ith.eq.1) then
        zwk(inthbdys,1:inrhos-1)=zwk(1,1:inrhos-1)
     endif
  enddo

  if(ier.eq.0) then
     do irho=1,inrhos-1
        do ith=1,inthbdys-1
           result(ith,irho)= &
                th_result(ith,irho+1) - th_result(ith,irho) +&
                zwk(ith+1,irho) - zwk(ith,irho)
        enddo
     enddo

  else
     result = 0
  endif

  deallocate(zwk,rr,zz,drdrho)

  call xpeval_free(rho_intrp)
  call xpeval_free(th_intrp)

end subroutine xplasma_integ_gwork

subroutine xplasma_integ_g10info(x5,w5)

  ! return normalized displacements and weights for 10 point Gauss
  ! non-adaptive integration evaluation

  implicit NONE

  real*8, intent(out) :: x5(5),w5(5)

  !----------------
  real*8 :: x1(5),w10(5)

  ! gauss-kronrod-patterson quadrature coefficients for use in
  ! quadpack routine qng.  these coefficients were calculated with
  ! 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.

  data x1    (  1) / 0.973906528517171720077964012084452D0/
  data x1    (  2) / 0.865063366688984510732096688423493D0/
  data x1    (  3) / 0.679409568299024406234327365114874D0/
  data x1    (  4) / 0.433395394129247190799265943165784D0/
  data x1    (  5) / 0.148874338981631210884826001129720D0/

  data w10   (  1) / 0.066671344308688137593568809893332D0/
  data w10   (  2) / 0.149451349150580593145776339657697D0/
  data w10   (  3) / 0.219086362515982043995534934228163D0/
  data w10   (  4) / 0.269266719309996355091226921569469D0/
  data w10   (  5) / 0.295524224714752870173892994651338D0/
  !----------------

  x5=x1
  w5=w10

end subroutine xplasma_integ_g10info

subroutine xplasma_integ_mkvec(nxbrk,xbrk,xvec,wx)

  ! generate evaluation vector for segmented integration using 10 point Gauss
  ! formula

  ! for integrating over splines, it is best to set the breaks at the
  ! spline node points, so that each sub interval integrated over is
  ! C-infinity...

  implicit NONE

  integer, intent(in) :: nxbrk        ! no. of x breaks (zone bdys)
  real*8, intent(in) :: xbrk(nxbrk) ! x break points

  real*8, intent(out) :: xvec(10*(nxbrk-1))  ! x evaluation points
  real*8, intent(out) :: wx(10*(nxbrk-1))    ! associated d(x) weights

  !---------------------------------------
  integer :: ixzon,ixp,ii
  real*8 :: zxcen,zdx2
  real*8, parameter :: HALF = 0.5d0
  !---------------------------------------

  real*8 :: x5(5),w5(5)

  call xplasma_integ_g10info(x5,w5)

  ixp = -10
  do ixzon = 1,nxbrk-1

     ixp=ixp+10

     zxcen= (xbrk(ixzon)+xbrk(ixzon+1))*HALF
     zdx2 = (xbrk(ixzon+1)-xbrk(ixzon))*HALF

     do ii=1,5
        xvec(ixp+ii)=zxcen+x5(ii)*zdx2
        xvec(ixp+5+ii)=zxcen-x5(ii)*zdx2
        wx(ixp+ii)=w5(ii)*zdx2
        wx(ixp+5+ii)=w5(ii)*zdx2
     enddo
  enddo

end subroutine xplasma_integ_mkvec

subroutine xplasma_integ_num(integ_name,integ_num,ilsurf, &
     dvol_need,dvdrho_need,bmaxa_need)

  use xplasma_definitions
  use eqi_intg_module
  implicit NONE

  character*(*), intent(in) :: integ_name
  integer, intent(out) :: integ_num
  logical, intent(out) :: ilsurf
  logical, intent(out) :: dvol_need
  logical, intent(out) :: dvdrho_need
  logical, intent(out) :: bmaxa_need

  !  translate integral/average character string name into integer id
  !  no errors-- return 0 if unrecognized.  CALLER should put name into
  !  UPPERCASE first!

  !---------------
  integer :: ilen,idtmp,iertmp
  !---------------

  integ_num=0

  ilen=max(2, len(trim(integ_name)) )

  !  common case: <1/R^2>S vs.  <1/R^2>...
  !    also support:  <fun> <fun>S I(fun) I(fun)S which respectively
  !    denote:  zonal volume average, 
  !             flux surface differential volume average,
  !             zonal integral, 
  !             surface integral
  !    to be applied to user specified function "fun" which will be checked
  !    for eligibility.

  bmaxa_need = .FALSE.

  if(integ_name(max(1,ilen-1):ilen).eq.'>S') then
     ! flux surface differential volume average
     ilsurf=.TRUE.
     ilen=ilen-1
     dvol_need=.FALSE.
     dvdrho_need=.TRUE.

  else if(integ_name(ilen:ilen).eq.'>') then
     ! zonal volume average
     ilsurf=.FALSE.
     dvol_need=.TRUE.
     dvdrho_need=.FALSE.

  else if(integ_name(max(1,ilen-1):ilen).eq.')S') then
     ! surface integral
     ilsurf=.TRUE.
     dvol_need=.FALSE.
     dvdrho_need=.FALSE.
     ilen=ilen-1

  else if(integ_name(ilen:ilen).eq.')') then
     ! zonal integral
     ilsurf=.FALSE.
     dvol_need=.FALSE.
     dvdrho_need=.FALSE.

  endif

  !  more specific names tested here...

  if(integ_name(1:ilen).eq.'VOL') then
     integ_num=eqi_intg_vol
     dvol_need=.TRUE.
     dvdrho_need=.FALSE.
     ilsurf=.TRUE.

  else if(integ_name(1:ilen).eq.'DVOL') then
     integ_num=eqi_intg_dvol
     dvol_need=.TRUE.
     dvdrho_need=.FALSE.
     ilsurf=.TRUE.  ! dVol will be differenced from VOL...

  else if(integ_name(1:ilen).eq.'AREA') then
     integ_num=eqi_intg_area
     dvdrho_need=.FALSE.
     dvol_need=.FALSE.
     ilsurf=.TRUE.

  else if(integ_name(1:ilen).eq.'DAREA') then
     integ_num=eqi_intg_darea
     dvdrho_need=.FALSE.
     dvol_need=.FALSE.
     ilsurf=.TRUE.  ! dAREA will be differenced from AREA...

  else if(integ_name(1:ilen).eq.'DVDRHO') then
     integ_num=eqi_intg_dvdrho
     dvdrho_need=.TRUE.
     dvol_need=.FALSE.
     ilsurf=.TRUE.

  else if(integ_name(1:ilen).eq.'ITOR')  then
     integ_num=eqi_intg_itor
     dvdrho_need=.FALSE.
     dvol_need=.FALSE.
     ilsurf=.TRUE.

  else if(integ_name(1:ilen).eq.'LPOL')  then
     integ_num=eqi_intg_lpol
     dvdrho_need=.FALSE.
     dvol_need=.FALSE.
     ilsurf=.TRUE.

  else if(integ_name(1:ilen).eq.'SURF')  then
     integ_num=eqi_intg_surf
     dvdrho_need=.FALSE.
     dvol_need=.FALSE.
     ilsurf=.TRUE.

  else if(integ_name(1:ilen).eq.'BOOZER_QOVG')  then
     ! integral related to Boozer "zeta"; loop integration -> q/g
     integ_num=eqi_intg_booz_qovg
     dvdrho_need=.FALSE.
     dvol_need=.FALSE.
     ilsurf=.TRUE.

  else if(integ_name(1:ilen).eq.'<1/R^2>') then
     integ_num=eqi_intg_r2i

  else if(integ_name(1:ilen).eq.'<GRAD(RHO)^2/R^2>') then
     integ_num=eqi_intg_g2r2i

  else if(integ_name(1:ilen).eq.'<1/R>') then
     integ_num=eqi_intg_ri

  else if(integ_name(1:ilen).eq.'<1/R^3>') then
     integ_num=eqi_intg_r3i

  else if(integ_name(1:ilen).eq.'<R^2>') then
     integ_num=eqi_intg_r2

  else if(integ_name(1:ilen).eq.'<R>') then
     integ_num=eqi_intg_r1

  else if(integ_name(1:ilen).eq.'<B*R^2>') then
     integ_num=eqi_intg_br2

  else if(integ_name(1:ilen).eq.'<1/B>') then
     integ_num=eqi_intg_bi

  else if(integ_name(1:ilen).eq.'<1/B^2>') then
     integ_num=eqi_intg_b2i

  else if(integ_name(1:ilen).eq.'<B^2>') then
     integ_num=eqi_intg_b2

  else if(integ_name(1:ilen).eq.'<B>') then
     integ_num=eqi_intg_b1

  else if(integ_name(1:ilen).eq.'<BZ^2>') then
     integ_num=eqi_intg_bz2

  else if(integ_name(1:ilen).eq.'<GRAD(RHO)^2/B^2>') then
     integ_num=eqi_intg_g2b2i

  else if(integ_name(1:ilen).eq.'<GRAD(RHO)>') then
     integ_num=eqi_intg_g1

  else if(integ_name(1:ilen).eq.'<GRAD(RHO)^2>') then
     integ_num=eqi_intg_g2

  else if(integ_name(1:ilen).eq.'<GRAD(RHO)^2/R^3>') then
     integ_num=eqi_intg_g2r3i

  else if(integ_name(1:ilen).eq.'<R^2*GRAD(RHO)^2>') then
     integ_num=eqi_intg_g2r2

  else if(integ_name(1:ilen).eq.'<1/(R*GRAD(RHO))>') then
     integ_num=eqi_intg_giri

  else if(integ_name(1:ilen).eq.'<F(H)NCLASS>') then
     bmaxa_need=.TRUE.
     integ_num = eqi_intg_fhnclass

  else if(integ_name(1:ilen).eq.'<F(H)NC_TSC>') then
     bmaxa_need=.TRUE.
     integ_num = eqi_intg_fhnc_tsc

  else if(integ_name(1:ilen).eq.'<SQRT(1-B/BMAX)>') then
     bmaxa_need=.TRUE.
     integ_num = eqi_intg_fhmmx

  else

     if((integ_name(1:1).eq.'<').and.(integ_name(ilen:ilen).eq.'>')) then
        idtmp=0
        if(ilen.gt.2) then
           !  see if item inside brackets is a known "integrable" profile
           !  for now, this means f(rho,theta) or f(theta,rho).

           call xplasma_integ_eligible(sp,integ_name(2:ilen-1),idtmp,iertmp)

        endif
        if(idtmp.eq.0) then
           call xplasma_errmsg_append(sp, &
                'cannot flux surface average '//trim(integ_name)// &
                ' -- not a profile.')
        else
           integ_num=eqi_intg_bign + idtmp
        endif

     else if((integ_name(1:2).eq.'I(').and.(integ_name(ilen:ilen).eq.')')) then
        idtmp=0
        if(ilen.gt.3) then
           !  see if item inside brackets is a known "integrable" profile
           !  for now, this means f(rho,theta) or f(theta,rho).

           call xplasma_integ_eligible(sp,integ_name(3:ilen-1),idtmp,iertmp)

        endif
        if(idtmp.eq.0) then
           call xplasma_errmsg_append(sp, &
                'cannot integrate '//trim(integ_name)// &
                ' -- not a profile.')
        else
           integ_num=2*eqi_intg_bign + idtmp
        endif

     endif

  endif

end subroutine xplasma_integ_num

subroutine xplasma_integ_mkindex(insurfs,rhos,infine,rhof,indexb)
  ! create index vector to point from fine grid to coarse grid
  ! (integration grid -> results grid)

  ! written here as if grid is "rho"; it could be "theta"...

  integer, intent(in) :: insurfs
  real*8, intent(in) :: rhos(insurfs)  ! coarse results grid

  integer, intent(in) :: infine
  real*8, intent(in) :: rhof(infine)   ! fine grid
  integer, intent(out) :: indexb(infine) ! indexing out

  !-----------------------
  integer :: i, cur_index
  !-----------------------

  cur_index=1

  do i=1,infine
     if(rhof(i).gt.rhos(cur_index+1)) cur_index = cur_index + 1

     indexb(i) = cur_index

  enddo
end subroutine xplasma_integ_mkindex
