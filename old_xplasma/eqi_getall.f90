subroutine eqi_getall(ivec,ccwflag, &
     zR,zZ,zphi,imap,zrho,zchi,tol,igrad, &
     ibflg,bvec,gbtensr, &
     nlist,inames,ivecd,fvals,fgrads, &
     nregion,ierr)

  !  **vectorized**
  !  get Bvector, grad(B) tensor, list of fcn values & list of fcn gradients
  !  at a particular spatial location

  use xplasma_obj_instance
  use eq_module

  IMPLICIT NONE

  integer ivec                        ! vector dimensioning, input
  logical :: ccwflag                  ! T for normal CCW poloidal angle
  ! F means transforms will be applied...

  REAL*8 zR(ivec),zZ(ivec),zphi(ivec) ! (R,Z,phi) point [*]

  integer imap                        ! control flag:
  !           =1 to map (R,Z,phi)->(rho,chi,phi)
  !            return vectors in (R,Z,phi) system, e.g. BR,BZ,Bphi
  !           =2 to map (rho,chi,phi)->(R,Z,phi)
  !            return vectors in (R,Z,phi) system, e.g. BR,BZ,Bphi
  !           =3 to map (rho,chi,phi)->(R,Z,phi)
  !            return vectors in (rho,chi,phi) system, e.g. Bperp=0, Bchi, Bphi

  !          see igrad (below) to determine if gradient is wanted.
  !          see ibflg (below) to determine if B is wanted

  !          see list of functions (below) to see what functions are wanted

  REAL*8 zrho(ivec),zchi(ivec)      ! (rho,chi,[phi]) magnetic flux coords
  !                                         [*]
  !  note that phi(...) is input in either case

  !  input:
  REAL*8 tol                        ! rel. tol. for (R,Z,phi)->flux coords

  !  input:
  integer igrad                     ! =1:  get gradients

  integer ibflg                     ! =1:  get B vector, grad(B) tensor
  !  grad(B) fetched only if igrad=ibflg=1

  !  output  (if ibflg=1):
  REAL*8 bvec(3,ivec)               ! (BR,BZ,Bphi) or (Bperp=0,Bpol,Bphi)
  REAL*8 gbtensr(3,3,ivec)          ! [grad(BR),grad(BZ),grad(Bphi)]

  !  here grad means ( d/dR, d/dZ, (1/R)*d/dphi )

  !  gbtensr(1:3,1)=grad(BR)
  !  gbtensr(1:3,2)=grad(BZ)
  !  gbtensr(1:3,3)=grad(Bphi)

  !  input:
  integer nlist                     ! number of functions to evaluate
  integer inames(*)                 ! function id numbers

  integer ivecd                     ! output vector dimensioning

  !  nlist=0 is allowed

  !  output (if nlist.gt.0):
  REAL*8 fvals(ivecd,*)             ! function values
  REAL*8 fgrads(3,ivecd,*)          ! function gradients

  !  output:
  integer nregion(ivec)             ! =1:  inside, .gt.1: outside core

  !  if input coordinates are (rho,chi,phi), nregion(j)=4 if rho is out of
  !  range of the core plasma, nregion(j)=1 otherwise; nregion(j)=5 for 
  !  rho(j) < 0.0

  !  if input coordinates are (R,Z,phi), nregion(i)=1 if (R,Z,phi)(i) is inside
  !  the core plasma; nregion(i)=2 if (R,Z,phi)(i) is outside the core but
  !  inside the limiters; nregion(i)=3 if (R,Z,phi)(i) is inside
  !  [Rmin,Rmax]x[Zmin,Zmax] but outside the limiters; nregion(i)=4 if
  !  (R,Z,phi)(i) is outside [Rmin,Rmax]x[Zmin,Zmax]

  integer ierr                      ! =0:  OK, .ne.0:  serious error.

  !  NOTE: field and profile evaluations are only for points (i) such that
  !  nregion(i).lt.4; if a profile only exists in f(rho,theta) form and
  !  imap.eq.1 and. nregion(i).gt.1, then, f(1,theta) is evaluated, where
  !  an extrapolated theta value has been chosen.  For points not evaluated,
  !  zeroes are returned for function and field values and gradients alike.
  !---------------------------------------------
  !  local:

  !  coordinates as used (forced to be in range)
  real*8 zrho_use(ivec),zR_use(ivec),zZ_use(ivec)
  real*8 zchi_use(ivec),zphi_use(ivec),zdist(ivec)

  integer :: ii,jj,kk,ifcn,jvec,imap_type

  real*8 :: zRmin,zRmax,zZmin,zZmax

  real*8, parameter :: ZERO=0.0d0
  real*8, parameter :: ONE =1.0d0

  integer id_fRZ(nlist+3),id_fmag(nlist+3),id_alt,igt

  integer :: neval_fRZ,neval_fmag
  integer :: ideval_fRZ(3*(nlist+3)),ideval_fmag(3*(nlist+3))
  integer :: indx,indx_fRZ(3*(nlist+3)),indx_fmag(3*(nlist+3))
  integer :: ideriv,ideriv_mag(3*(nlist+3),2),ideriv_RZ(3*(nlist+3),2)

  real*8, dimension(:), allocatable :: drhodR,drhodZ,dthdR,dthdZ
  real*8, dimension(:,:), allocatable :: zwk
  real*8 :: elrho(2),elth(2),zdenom,zperp,zpll,zdfdr,zdfdz,bdytol
  real*8 :: zBR,zBZ,zBpol,zBpolgR,zBpolgZ

  !---------------------------------------------
  !  eqi_getall "old_xplasma" conversion dmc Jul 2006

  integer :: id_Bmod,id_BR,id_BZ,id_R,id_Z,id_g,itype_lim,itype_obj,iertmp
  integer :: nsnccwB,nsnccwI,nsnccw_th
  logical :: axisymm
  character*32 zname

  logical :: RZpref,RZout,gtrans

  !---------------------------------------------

  RZpref = (imap.eq.1)  ! prefer f(R,Z) evaluation if (R,Z) are input
  !                       otherwise prefer f(rho,theta) evaluation
  !  also means (R,Z) are input, not (rho,theta)

  RZout = (imap.le.2)   ! return [R,Z] gradients & BR, BZ
  !                       otherwise: [perp,pll] gradients & Btheta

  !.... set default outputs

  if(RZpref) then
     zrho=ZERO
     zchi=ZERO
  else
     zR=ZERO
     zZ=ZERO
  endif
  nregion=0

  if((igrad.gt.1).or.(igrad.lt.0)) then
     write(lunerr,*) ' ?eqi_getall: igrad=',igrad
     write(lunerr,*) '  only a value of 0 or 1 should be given.'
     ierr=9999
  endif

  if((ibflg.gt.1).or.(ibflg.lt.0)) then
     write(lunerr,*) ' ?eqi_getall: ibflg=',ibflg
     write(lunerr,*) '  only a value of 0 or 1 should be given.'
     ierr=9999
  endif

  if(ibflg.eq.1) then
     bvec=ZERO
     if(igrad.eq.1) gbtensr=ZERO
  endif

  ierr=0
  if(nlist.lt.0) then
     write(lunerr,*) ' ?eqi_getall: nlist=',nlist
     write(lunerr,*) '  a non-negative value is required.'
     ierr=9999
  else if(nlist.gt.0) then
     fvals(1:ivec,1:nlist)=ZERO
     if(igrad.eq.1) fgrads(1:3,1:ivec,1:nlist)=ZERO
  endif
  if(ierr.ne.0) return

  !-----------------------------------------------------------
  !  now do error checking (other than nlist.ge.0):

  if(ccwflag) then
     nsnccw_th=-1   ! CCW j -> Bpol direction opposite to d[R,Z]/dtheta
  else
     nsnccw_th=1
  endif

  call xplasma_global_info(s,ierr,axisymm=axisymm, &
       bphi_ccw=nsnccwB,jphi_ccw=nsnccwI,bdytol=bdytol)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqi_getall: unexpected xplasma_global_info error:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif
  if(.not.axisymm) then
     call eq_errmsg(' eqi_getall:  non-axisymmetry not supported!')
     ierr=1
     return
  endif

  call xplasma_common_ids(s,ierr, id_Bmod=id_Bmod,id_BR=id_BR,id_BZ=id_BZ, &
       id_g=id_g,id_R=id_R,id_Z=id_Z)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqi_getall: unexpected xplasma_common_ids error:'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  if(min(id_R,id_Z).eq.0) then
     ierr=9999
     write(lunerr,*) ' ?eqi_getall: equilibrium geometry not ready.'
  endif

  if(min(id_Bmod,id_BR,id_BZ).eq.0) then
     ierr=9999
     write(lunerr,*) ' ?eqi_getall: equilibrium B field not ready.'
  endif

  if(ierr.ne.0) return

  if(imap.eq.1) then
     call xplasma_lim_info(s,ierr,itype=itype_lim)
     if(ierr.ne.0) then
        itype_lim=0
        call xplasma_error(s,ierr,lunerr)
     endif
     if(itype_lim.eq.0) then
        write(lunerr,*) ' ?eqi_getall: limiter specification not found.'
        ierr=9999
     endif
  endif

  do ifcn=1,nlist
     id_fRZ(ifcn)=0
     id_fmag(ifcn)=0

     call xplasma_get_item_info(s,inames(ifcn),iertmp,itype=itype_obj, &
          name=zname)

     if((iertmp.ne.0).or.(itype_obj.eq.0)) then
        if(iertmp.ne.0) call xplasma_error(s,iertmp,lunerr)
        write(lunerr,*) ' ?eqi_getall: invalid xplasma id code:'
        write(lunerr,*) '  inames(',ifcn,')=',inames(ifcn)
        ierr=9999

     else if(itype_obj.ne.xplasma_profType) then
        write(lunerr,*) ' ?eqi_getall: not a profile id code:'
        write(lunerr,*) '  inames(',ifcn,')=',inames(ifcn)
        write(lunerr,*) '    name       = "',trim(zname),'".'
        ierr=9999

     else
        ! OK if grid type is recognized...

        call xplasma_prof_info(s,inames(ifcn),iertmp, &
             gridType=igt, profId1=id_alt)

        if(igt.eq.1) then
           id_fmag(ifcn)=inames(ifcn)
           if(id_alt.gt.0) id_fRZ(ifcn)=id_alt

        else if(igt.eq.2) then
           id_fRZ(ifcn)=inames(ifcn)
           if(id_alt.gt.0) id_fmag(ifcn)=id_alt

        else
           ierr=9999
           write(lunerr,*) ' ?eqi_getall: profile "',trim(zname),'" is not'
           write(lunerr,*) '  defined over a supported set of axes;'
           write(lunerr,*) '  only functions of form f(R,Z) and f(rho,theta) are allowed.'
        endif
     endif

  enddo

  if(ierr.ne.0) return

  !------------------------------------------------

  if(ibflg.eq.1) then
     id_fmag(nlist+1)=id_BR
     id_fRZ(nlist+1)=id_BR  ! eval routine will convert
     id_fmag(nlist+2)=id_BZ
     id_fRZ(nlist+2)=id_BZ  ! eval routine will convert

     id_fmag(nlist+3)=id_g
     id_fRZ(nlist+3)=0      ! use g(rho)/R
  else
     id_fmag(nlist+1:nlist+3)=0
     id_fRZ(nlist+1:nlist+3)=0
  endif

  !------------------------------------------------
  !  coordinate mapping...

  call xplasma_RZminmax_extended(s, zRmin,zRmax, zZmin,zZmax, ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqi_getall: unexpected xplasma_RZminmax_extended error!'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  !----------------------------------
  if(RZpref) then

     !  compute (R,Z,phi) -> (rho,theta,phi)

     if(tol.ge.1.0d-3) then
        imap_type=xp_imap_polar              ! approximation OK
     else
        imap_type=xp_imap_newton             ! use slower Newton iteration
     endif

     call xplasma_ctrans(s,xp_use_coord_vecsize,ierr, &
          r_in=zR, z_in=zZ, tol=tol, maptype=imap_type, &
          rho_out=zrho, theta_out=zchi, nregion=nregion, ccw_theta=ccwflag)
     if(ierr.ne.0) then
        write(lunerr,*) ' ?eqi_getall: unexpected xplasma_ctrans(R,Z) error!'
        call xplasma_error(s,ierr,lunerr)
        return
     endif

     !  for points outside, set up nregion category appropriately

     jvec=0
     do ii=1,ivec
        if(nregion(ii).eq.3) then
           nregion(ii)=4
        else if((zR(ii).lt.zRmin).or.(zR(ii).gt.zRmax).or. &
             (zZ(ii).lt.zZmin).or.(zZ(ii).gt.zZmax)) then
           nregion(ii)=4
        endif

        if(nregion(ii).eq.2) then
           jvec=jvec+1
           zR_use(jvec)=zR(ii)
           zZ_use(jvec)=zZ(ii)
        endif
     enddo

     if(jvec.gt.0) then
        call xplasma_lim_distance(s,zR_use(1:jvec),zZ_use(1:jvec), &
             zdist(1:jvec),ierr, maptype=imap_type)
        if(ierr.ne.0) then
           write(lunerr,*) ' ?eqi_getall: unexpected xplasma_lim_distance error!'
           call xplasma_error(s,ierr,lunerr)
           return
        endif
     endif

     jvec=0
     do ii=1,ivec
        if(nregion(ii).eq.2) then
           jvec=jvec+1
           if(zdist(jvec).ge.0.0d0) then
              nregion(ii)=3
           endif
        endif
     enddo

  !----------------------------------
  else

     !  compute (rho,theta,phi) -> (R,Z,phi)

     jvec=0
     do ii=1,ivec
        if((zrho(ii).gt.-bdytol).and.(zrho(ii).lt.ONE+bdytol)) then
           nregion(ii)=1
           jvec=jvec+1
           zrho_use(jvec)=max(ZERO,min(ONE,zrho(ii)))
           zchi_use(jvec)=zchi(ii)

        else if(zrho(ii).le.-bdytol) then
           nregion(ii)=5
           zR(ii)=ZERO
           zZ(ii)=ZERO

        else if(zrho(ii).ge.ONE+bdytol) then
           nregion(ii)=4
           zR(ii)=ZERO
           zZ(ii)=ZERO

        endif
     enddo

     if(jvec.gt.0) then
        call xplasma_ctrans(s,xp_use_coord_vecsize,ierr, &
             rho_in=zrho_use(1:jvec),theta_in=zchi_use(1:jvec), &
             R_out=zR_use(1:jvec),Z_out=zZ_use(1:jvec), ccw_theta=ccwflag)
        if(ierr.ne.0) then
           write(lunerr,*) ' ?eqi_getall: unexpected xplasma_ctrans(rho,theta) error!'
           call xplasma_error(s,ierr,lunerr)
           return
        endif
     endif

     jvec=0
     do ii=1,ivec
        if(nregion(ii).eq.1) then
           jvec=jvec+1
           zR(ii)=zR_use(jvec)
           zZ(ii)=zZ_use(jvec)
        endif
     enddo
  endif

  !-------------------------------------------
  !  function and field evaluation...
  !  for each (id) choose f(R,Z) or f(rho,theta) -- the one available, or,
  !  the one preferred if there is a choice (prior error checking assures
  !  that at least one is available).
        
  neval_fRZ=0
  neval_fmag=0

  do ifcn=1,nlist+3
     if(RZpref) then
        if(id_fRZ(ifcn).gt.0) then
           call add_fRZ
        else
           call add_fmag
        endif
     else
        if(id_fmag(ifcn).gt.0) then
           call add_fmag
        else
           call add_fRZ
        endif
     endif
  enddo

  !----------------------
  ! evaluation points

  jvec=0
  do ii=1,ivec
     if(nregion(ii).lt.4) then
        jvec=jvec+1
        zR_use(jvec)=zR(ii)
        zZ_use(jvec)=zZ(ii)
        zrho_use(jvec)=max(ZERO,min(ONE,zrho(ii)))
        zchi_use(jvec)=zchi(ii)
     endif
  enddo

  if(jvec.eq.0) return

  if(igrad.eq.0) then
     gtrans = .FALSE.
  else
     gtrans = (neval_fmag.gt.0).or.(.not.RZout)
  endif

  ierr=0
  if(gtrans) then
     allocate(drhodR(jvec),drhodZ(jvec),dthdR(jvec),dthdZ(jvec))
     call xplasma_rzjac(s,zrho_use(1:jvec),zchi_use(1:jvec),ierr, &
          drhodr=drhodR,drhodz=drhodZ, dthdr=dthdR,dthdz=dthdZ, &
          ccwflag=ccwflag)
     if(ierr.ne.0) then
        deallocate(drhodR,drhodZ,dthdR,dthdZ)
        return
     endif
  endif
  
  if(neval_fRZ.gt.0) then
     allocate(zwk(jvec,neval_fRZ))
     call xplasma_eval_prof(s,ideval_fRZ(1:neval_fRZ), &
          xplasma_R_coord,zR_use(1:jvec),xplasma_Z_coord,zZ_use(1:jvec), &
          zwk, iertmp, &
          ideriv1s=ideriv_RZ(1:neval_fRZ,1), &
          ideriv2s=ideriv_RZ(1:neval_fRZ,2))
     if(iertmp.ne.0) then
        ierr=iertmp
        write(lunerr,*) ' ?eqi_getall: unexpected xplasma_eval(R,Z) error!'
        call xplasma_error(s,ierr,lunerr)
     else
        jvec=0
        do ii=1,ivec
           if(nregion(ii).lt.4) then
              jvec=jvec+1

              !  indx.eq.-3 will not arise because Bphi is evaluated
              !  using nsnccwB*g(rho)/R (not f(R,Z))

              do jj=1,neval_fRZ
                 indx=indx_fRZ(jj)
                 if(igrad.eq.0) then
                    if(indx.gt.0) then
                       fvals(ii,indx)=zwk(jvec,jj)
                    else if(indx.eq.-1) then
                       bvec(1,ii)=zwk(jvec,jj)
                    else if(indx.eq.-2) then
                       bvec(2,ii)=zwk(jvec,jj)
                    endif
                 else
                    ideriv=1*ideriv_RZ(jj,1) + 2*ideriv_RZ(jj,2)
                    if(indx.gt.0) then
                       if(ideriv.eq.0) then
                          fvals(ii,indx)=zwk(jvec,jj)
                       else if(ideriv.eq.1) then
                          fgrads(1,ii,indx)=zwk(jvec,jj)
                       else if(ideriv.eq.2) then
                          fgrads(2,ii,indx)=zwk(jvec,jj)
                       endif
                    else if(indx.eq.-1) then
                       if(ideriv.eq.0) then
                          bvec(1,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.1) then
                          gbtensr(1,1,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.2) then
                          gbtensr(2,1,ii)=zwk(jvec,jj)
                       endif
                    else if(indx.eq.-2) then
                       if(ideriv.eq.0) then
                          bvec(2,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.1) then
                          gbtensr(1,2,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.2) then
                          gbtensr(2,2,ii)=zwk(jvec,jj)
                       endif
                    endif
                 endif
              enddo
              if((igrad.eq.1).and.(.not.RZout)) then
                 ! transform gradients --> [perp,pll]
                 zdenom=sqrt(drhodR(jvec)**2+drhodZ(jvec)**2)
                 elrho(1)=drhodR(jvec)/zdenom
                 elrho(2)=drhodZ(jvec)/zdenom
                 if(ccwflag) then
                    elth(1)=-elrho(2); elth(2)= elrho(1)
                 else
                    elth(1)= elrho(2); elth(2)=-elrho(1)
                 endif
                 do jj=1,neval_fRZ,3
                    indx=indx_fRZ(jj)
                    if(indx.gt.0) then
                       zperp=elrho(1)*fgrads(1,ii,indx) + &
                             elrho(2)*fgrads(2,ii,indx)
                       zpll =elth(1)*fgrads(1,ii,indx) + &
                             elth(2)*fgrads(2,ii,indx)
                       fgrads(1,ii,indx)=zperp
                       fgrads(2,ii,indx)=zpll
                    endif
                 enddo
              endif
           endif
        enddo
     endif
     deallocate(zwk)
  endif

  if(neval_fmag.gt.0) then
     allocate(zwk(jvec,neval_fmag))
     call xplasma_eval_prof(s,ideval_fmag(1:neval_fmag), &
          xplasma_rho_coord,zrho_use(1:jvec), &
          xplasma_theta_coord,zchi_use(1:jvec), &
          zwk, iertmp, ccwflag2=ccwflag, &
          ideriv1s=ideriv_mag(1:neval_fmag,1), &
          ideriv2s=ideriv_mag(1:neval_fmag,2))
     if(iertmp.ne.0) then
        ierr=iertmp
        write(lunerr,*) ' ?eqi_getall: unexpected xplasma_eval(rho,theta) error!'
        call xplasma_error(s,ierr,lunerr)
     else
        jvec=0
        do ii=1,ivec
           if(nregion(ii).lt.4) then
              jvec=jvec+1

              !  evaluation gives {f,df/drho,df/dtheta}
              !  gradients to be built later--
              !    except Bphi = nsnccwB*g/R
              !    and grad(Bphi) = nsnccwB*(dg/drho)*grad(rho)/R - g/R**2
              !      are evaluated directly

              do jj=1,neval_fmag
                 indx=indx_fmag(jj)
                 if(igrad.eq.0) then
                    if(indx.gt.0) then
                       fvals(ii,indx)=zwk(jvec,jj)
                    else if(indx.eq.-1) then
                       bvec(1,ii)=zwk(jvec,jj)
                    else if(indx.eq.-2) then
                       bvec(2,ii)=zwk(jvec,jj)
                    else if(indx.eq.-3) then
                       bvec(3,ii)=nsnccwB*zwk(jvec,jj)/zR_use(jvec)
                    endif
                 else
                    ideriv=1*ideriv_mag(jj,1) + 2*ideriv_mag(jj,2)
                    if(indx.gt.0) then
                       if(ideriv.eq.0) then
                          fvals(ii,indx)=zwk(jvec,jj)
                       else if(ideriv.eq.1) then
                          fgrads(1,ii,indx)=zwk(jvec,jj)
                       else if(ideriv.eq.2) then
                          fgrads(2,ii,indx)=zwk(jvec,jj)
                       endif
                    else if(indx.eq.-1) then
                       !  BR
                       if(ideriv.eq.0) then
                          bvec(1,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.1) then
                          gbtensr(1,1,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.2) then
                          gbtensr(2,1,ii)=zwk(jvec,jj)
                       endif
                    else if(indx.eq.-2) then
                       !  BZ
                       if(ideriv.eq.0) then
                          bvec(2,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.1) then
                          gbtensr(1,2,ii)=zwk(jvec,jj)
                       else if(ideriv.eq.2) then
                          gbtensr(2,2,ii)=zwk(jvec,jj)
                       endif
                    else if(indx.eq.-3) then
                       !  Bphi
                       if(ideriv.eq.0) then
                          bvec(3,ii)=nsnccwB*zwk(jvec,jj)/zR_use(jvec)
                       else if(ideriv.eq.1) then
                          gbtensr(1,3,ii)=nsnccwB*zwk(jvec,jj)*drhodR(jvec) -&
                               bvec(3,ii)/zR_use(jvec)
                          gbtensr(2,3,ii)=nsnccwB*zwk(jvec,jj)*drhodZ(jvec)
                       endif
                       !  (ideriv.eq.2 refers to theta derivative = 0 here)
                    endif
                 endif
              enddo
              if(igrad.eq.1) then

                 ! convert {df/drho,df/dtheta} to grad(f)

                 do jj=1,neval_fmag,3
                    indx=indx_fmag(jj)
                    if(indx.gt.0) then
                       zdfdr=fgrads(1,ii,indx)*drhodr(jvec) + &
                             fgrads(2,ii,indx)*dthdr(jvec)
                       zdfdz=fgrads(1,ii,indx)*drhodz(jvec) + &
                             fgrads(2,ii,indx)*dthdz(jvec)
                       fgrads(1,ii,indx)=zdfdr
                       fgrads(2,ii,indx)=zdfdz
                    else if((indx.eq.-1).or.(indx.eq.-2)) then
                       kk=-indx
                       zdfdr=gbtensr(1,kk,ii)*drhodr(jvec) + &
                             gbtensr(2,kk,ii)*dthdr(jvec)
                       zdfdz=gbtensr(1,kk,ii)*drhodz(jvec) + &
                             gbtensr(2,kk,ii)*dthdz(jvec)
                       gbtensr(1,kk,ii)=zdfdr
                       gbtensr(2,kk,ii)=zdfdz
                    endif
                 enddo

                 if(.not.RZout) then
                    ! transform gradients --> [perp,pll]
                    zdenom=sqrt(drhodR(jvec)**2+drhodZ(jvec)**2)
                    elrho(1)=drhodR(jvec)/zdenom
                    elrho(2)=drhodZ(jvec)/zdenom
                    if(ccwflag) then
                       elth(1)=-elrho(2); elth(2)= elrho(1)
                    else
                       elth(1)= elrho(2); elth(2)=-elrho(1)
                    endif
                    do jj=1,neval_fmag,3
                       indx=indx_fmag(jj)
                       if(indx.gt.0) then
                          zperp=elrho(1)*fgrads(1,ii,indx) + &
                               elrho(2)*fgrads(2,ii,indx)
                          zpll =elth(1)*fgrads(1,ii,indx) + &
                               elth(2)*fgrads(2,ii,indx)
                          fgrads(1,ii,indx)=zperp
                          fgrads(2,ii,indx)=zpll
                       endif
                    enddo
                 endif
              endif
           endif
        enddo
     endif
     deallocate(zwk)
  endif

  if((.not.RZout).and.(ierr.eq.0).and.(ibflg.eq.1)) then
     jvec=0
     do ii=1,ivec
        if(nregion(ii).lt.4) then
           jvec=jvec+1
           zBR=bvec(1,ii)
           zBZ=bvec(2,ii)
           zBpol=nsnccw_th*nsnccwI*sqrt(zBR**2+zBZ**2)
           bvec(1,ii)=ZERO
           bvec(2,ii)=zBpol
           if(igrad.eq.1) then
              zdenom=sqrt(drhodR(jvec)**2+drhodZ(jvec)**2)
              elrho(1)=drhodR(jvec)/zdenom
              elrho(2)=drhodZ(jvec)/zdenom
              if(ccwflag) then
                 elth(1)=-elrho(2); elth(2)= elrho(1)
              else
                 elth(1)= elrho(2); elth(2)=-elrho(1)
              endif

              !  grad(Bphi) in perp,pll frame...
              zperp=elrho(1)*gbtensr(1,3,ii) + &
                   elrho(2)*gbtensr(2,3,ii)
              zpll =elth(1)*gbtensr(1,3,ii) + &
                   elth(2)*gbtensr(2,3,ii)
              gbtensr(1,3,ii)=zperp
              gbtensr(2,3,ii)=zpll

              !  (R,Z) grad(BR) & grad(BZ) --> (R,Z) grad(Bpol)

              zBpolgR = (zBR*gbtensr(1,1,ii)+zBZ*gbtensr(1,2,ii))/zBpol
              zBpolgZ = (zBR*gbtensr(2,1,ii)+zBZ*gbtensr(2,2,ii))/zBpol

              !  convert to perp,pll frame

              zperp = elrho(1)*zBpolgR + elrho(2)*zBpolgZ
              zpll  = elth(1)*zBpolgR + elth(2)*zBpolgZ
              gbtensr(1,1,ii)=ZERO
              gbtensr(2,1,ii)=ZERO
              gbtensr(1,2,ii)=zperp
              gbtensr(2,2,ii)=zpll
           endif
        endif
     enddo
  endif

  if(allocated(drhodR)) deallocate(drhodR,drhodZ,dthdR,dthdZ)

  contains
    subroutine add_fRZ

      integer :: n

      if(igrad.eq.0) then
         neval_fRZ=neval_fRZ+1
         ideval_fRZ(neval_fRZ)=id_fRZ(ifcn)
         if(ifcn.le.nlist) then
            indx_fRZ(neval_fRZ)=ifcn
         else
            indx_fRZ(neval_fRZ)=-(ifcn-nlist)
         endif
         ideriv_RZ(neval_fRZ,1:2)=0
      else
         n=neval_fRZ
         ideval_fRZ(n+1:n+3)=id_fRZ(ifcn)
         if(ifcn.le.nlist) then
            indx_fRZ(n+1:n+3)=ifcn
         else
            indx_fRZ(n+1:n+3)=-(ifcn-nlist)
         endif
         ideriv_RZ(n+1:n+3,1:2)=0
         ideriv_RZ(n+2,1)=1
         ideriv_RZ(n+3,2)=1
         neval_fRZ=n+3
      endif

    end subroutine add_fRZ

    subroutine add_fmag

      integer :: n

      if(igrad.eq.0) then
         neval_fmag=neval_fmag+1
         ideval_fmag(neval_fmag)=id_fmag(ifcn)
         if(ifcn.le.nlist) then
            indx_fmag(neval_fmag)=ifcn
         else
            indx_fmag(neval_fmag)=-(ifcn-nlist)
         endif
         ideriv_mag(neval_fmag,1:2)=0
      else
         n=neval_fmag
         ideval_fmag(n+1:n+3)=id_fmag(ifcn)
         if(ifcn.le.nlist) then
            indx_fmag(n+1:n+3)=ifcn
         else
            indx_fmag(n+1:n+3)=-(ifcn-nlist)
         endif
         ideriv_mag(n+1:n+3,1:2)=0
         ideriv_mag(n+2,1)=1
         ideriv_mag(n+3,2)=1
         neval_fmag=n+3
      endif

    end subroutine add_fmag

end subroutine eqi_getall
