subroutine eq_flxint_arr2d(iorder,irhoax,zdata,id1,id2, &
     ibcrho0,idz0,zbcrho0,ibcrho1,idz1,zbcrho1, &
     iwant,result,inumchi,inumrho,ierr)

  ! integration of user defined data -- legacy f77 interface

  use xplasma_obj_instance
  use eq_module
  IMPLICIT NONE

  !  this routine will numerically integrate a 2d axisymmetric function
  !  based on array data provided by the user.  This is a mulltistep process:

  !     (1) a spline (or Hermite) representation of the array data is created
  !         (eqm_frhochi call)
  !     (2) this representation is integrated
  !         (fluxav_zonint or fluxav_zonavg or fluxav_surfint call)
  !     (3) the spline (or Hermite) representation is deleted.

  !         the argument iorder specifies the representation:
  !               iorder=2 --> spline
  !               iorder=1 --> Akima Hermite
  !               iorder=0 --> piecewise linear

  !  currently only axisymmetric data is supported, so the data is a 2d
  !  array, either:
  !     vs. (rho,chi), set irhoax=1; presume zdata(1:nrho,1:nchi), or,
  !     vs. (chi,rho), set irhoax=2; presume zdata(1:nchi,1:nrho)

  !  (rho is the radial flux coordinate; chi is the poloidal angle coordinate)
  !  these are the grids defined by prior eqm_rho and eqm_chi calls.

  !  note the actual size of first dimension of the data is "id1", which must
  !  be greater than either nrho or nchi, depending on the value of irhoax.

  !  ibcrho0,zbcrho0,ibcrho1,zbcrho1 -- boundary conditions for spline
  !  or Hermite fit (passed to eqm_frhochi)

  !  *** iwant *** specifies the type of integration desired:
  !    here "f" denotes the interpolating function generated from the
  !    passed array data.

  !    iwant=1:  surface integral:
  !       int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * f }

  !    iwant=5:  volume weighted surface integral:
  !        int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * f * det(Jacobian) }

  !      **caution** if rho(1) corresponds to the mag. axis, a
  !        det(Jacobian)=0 and a value of 0.0 is always returned!

  !    iwant=6:  volume weighted surface average:

  !        int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * f * det(Jacobian)}
  !        -----------------------------------------------
  !        int[phi=0 to 2pi] dphi * {
  !            int[chi=0 to 2pi] dchi * det(Jacobian)}

  !      **caution** if rho(1) corresponds to the mag. axis,
  !        the value "f" evaluated at the axis is returned there.

  !  zonal integrals:  result(1:nrho-1) are computed; result(j)
  !                    corresponds to the integral from rho(j) to rho(j+1)

  !    iwant=2:  integrate f*drho*dtheta*dphi:

  !       int[rho=rho(j) to rho(j+1)] drho * {
  !          int[phi=0 to 2pi] dphi * {
  !               int[chi=0 to 2pi] dchi * f } }

  !    iwant=3:  integrate f*dV (use det(Jacobian))

  !       int[rho=rho(j) to rho(j+1)] drho * {
  !          int[phi=0 to 2pi] dphi * {
  !               int[chi=0 to 2pi] dchi * f * det(Jacobian) } }

  !    iwant=4:  compute <f> = (integrated f*dV)/dVol
  !              where dVol is the volume btw rho(j) and rho(j+1)

  !  **caution** in axisymmetric case there is no variation with phi,
  !  but a factor of 2pi still arises from the phi integration...

  !----------------------------------
  !  arguments:
  !  -- input --

  integer iorder                ! order of fit (2=spline,1=Hermite,...)
  integer irhoax                ! which axis is "rho" (1 or 2)
  integer id1                   ! size of 1st dimension of zdata
  integer id2                   ! size of 2nd dimension of zdata
  real*8 zdata(id1,id2)         ! the data itself...

  !  irhoax=1 --> expect id1.ge.nrho, id2.ge.nchi
  !  irhoax=2 --> expect id1.ge.nchi, id2.ge.nrho

  integer idz0,idz1             ! sizes of BC data arrays

  integer ibcrho0,ibcrho1       ! BC type @rho(axis), @rho(bdy) resp.
  real*8 zbcrho0(idz0),zbcrho1(idz1)  ! corresponding BC data

  !  BC type values:  for ibcrho0:
  !     =0 -- use "not a knot", zbcrho0(...) ignored
  !     =1 -- match slope, given at x(1),chi(ichi) by zbcrho0(ichi)
  !     =2 -- match 2nd deriv., given at x(1),chi(ichi) by zbcrho0(ichi)
  !     =3 -- boundary condition is slope=0 (df/dx=0) at x(1), all chi(j)
  !           (zbcrho0(...) ignored).
  !  and similarly for (ibcrho1,zbcrho1(...)).

  integer iwant                 ! type of integration desired
  !                                   ! see description, above...

  integer inumchi              ! 1st (chi) dim of result
  integer inumrho              ! 2nd (rho) dim of result

  !  if inumchi=1, do integrations over the entire flux zone / flux surface
  !  if inumchi>1, do integrations over chi segments.

  !  -- output --

  real*8 result(inumchi,inumrho) ! results of integrations...

  integer ierr                  ! completion code, 0=normal

  !---------------------------------------------------------------------

  integer id_ax1,id_ax2,id_tmp,idi

  integer :: nchi,nrho,iertmp,inrhoi,inthi,inthzons,iadj
  integer :: id_R,id_rho,id_chi

  real*8, dimension(:), allocatable :: dvol_1d,dvdrho_1d
  real*8, dimension(:,:), allocatable :: dvol_2d,dvdrho_2d
  integer :: iwarn_vol,iwarn_dv

  character*10 zname
  logical :: ineed_dvol,ineed_dvdrho

  !---------------------------------------------------------------------
  !  quick checks...

  ierr=0
  result=0

  call xplasma_common_ids(s,ierr, id_R=id_R)
  if(ierr.ne.0) then
     return
  else if(id_R.eq.0) then
     return
  endif

  call xplasma_prof_info(s,id_R,iertmp, gridId1=id_chi, gridId2=id_rho)
  call xplasma_grid_size(s,id_chi,nchi,iertmp)
  call xplasma_grid_size(s,id_rho,nrho,iertmp)

  if((irhoax.lt.1).or.(irhoax.gt.2)) then
     write(lunerr,*) ' ?eq_flxint_arr2d: invalid: irhoax=',irhoax
     ierr=1
  endif
  if(ierr.ne.0) return

  if(irhoax.eq.1) then
     if(id1.lt.nrho) then
        write(lunerr,*) ' ?eq_flxint_arr2d:  id1=',id1,' too small.'
        write(lunerr,*) '  was expecting id1.ge.nrho, nrho = ',nrho
        ierr=1
     endif
     if(id2.lt.nchi) then
        write(lunerr,*) ' ?eq_flxint_arr2d:  id2=',id2,' too small.'
        write(lunerr,*) '  was expecting id2.ge.nchi, nchi = ',nchi
        ierr=1
     endif

  else
     if(id1.lt.nchi) then
        write(lunerr,*) ' ?eq_flxint_arr2d:  id1=',id1,' too small.'
        write(lunerr,*) '  was expecting id1.ge.nchi, nchi = ',nchi
        ierr=1
     endif
     if(id2.lt.nrho) then
        write(lunerr,*) ' ?eq_flxint_arr2d:  id2=',id2,' too small.'
        write(lunerr,*) '  was expecting id2.ge.nrho, nrho = ',nrho
        ierr=1
     endif
  endif

  if((ibcrho0.eq.1).or.(ibcrho0.eq.2)) then
     if(idz0.lt.nchi) then
        write(lunerr,*) ' ?eq_flxint_arr2d: idz0=',idz0,' too small'
        write(lunerr,*) '  was expecting idz0.ge.nchi, nchi = ',nchi
        write(lunerr,*) '  because ibcrho0=',ibcrho0
        ierr=1
     endif
  endif

  if((ibcrho1.eq.1).or.(ibcrho1.eq.2)) then
     if(idz1.lt.nchi) then
        write(lunerr,*) ' ?eq_flxint_arr2d: idz1=',idz1,' too small'
        write(lunerr,*) '  was expecting idz1.ge.nchi, nchi = ',nchi
        write(lunerr,*) '  because ibcrho1=',ibcrho1
        ierr=1
     endif
  endif

  if(ierr.ne.0) return

  !------------------------------------------
  !  (1)  set up interpolating function...

  if(irhoax.eq.1) then
     id_ax1=id_rho
     id_ax2=id_chi
  else
     id_ax1=id_chi
     id_ax2=id_rho
  endif

  call eqm_frhochi(iorder,id_ax1,id_ax2,'__TMPI',zdata,id1, &
       ibcrho0,zbcrho0,ibcrho1,zbcrho1,id_tmp,ierr)
  if(ierr.ne.0) return

  !  (2)  perform the integration

  ineed_dvol=.FALSE.
  ineed_dvdrho=.FALSE.

  if(iwant.eq.1) then

     iadj=0 ! surface oriented
     zname = 'I(__TMPI)S'

  else if(iwant.eq.2) then

     !  integral of f*dtheta*dphi*drho over zones

     iadj=1 ! zone oriented
     zname = 'I(__TMPI)'

  else if(iwant.eq.3) then

     !  integral of f*dtheta*dphi*drho*det(J) over zones

     iadj=1 ! zone oriented
     zname = '<__TMPI>'
     ineed_dvol = .TRUE.

  else if(iwant.eq.4) then

     !  integral of f*dtheta*dphi*drho*det(J) over zones, divided by
     !  integral of   dtheta*dphi*drho*det(J) over the same zones.

     iadj=1 ! zone oriented
     zname = '<__TMPI>'

  else if(iwant.eq.5) then

     !  integral of f*dtheta*dphi*det(J) on surfaces

     iadj=0 ! surface oriented
     zname = '<__TMPI>S'
     ineed_dvdrho = .TRUE.

  else if(iwant.eq.6) then

     !  integral of f*dtheta*dphi*det(J) on surfaces divided by
     !    integral of dtheta*dphi*det(J) on surfaces; on axis just return f.

     iadj=0 ! surface oriented
     zname = '<__TMPI>S'

  endif

  ! get 1d or 2d integrator dataset information:

  call eq_find_intg(inumchi,idi,inrhoi,inthi,ierr)
  if(ierr.ne.0) return

  if(inumchi.eq.1) then
     !  integrate over entire flux surface or zone

     if(inumrho.lt.inrhoi-iadj) then
        write(lunerr,*) ' ?eq_flxint: "inumrho" dimension too small.'
        write(lunerr,*) '  found: ',inumrho,' need at least: ',inrhoi-iadj
        ierr=1
        return
     endif
        
     allocate(dvol_1d(inrhoi-1),dvdrho_1d(inrhoi))

     call xplasma_rho_zonint(s,idi,zname,result(1,1:inrhoi-iadj),ierr, &
          dvol_out=dvol_1d,iwarn_dvol=iwarn_vol, &
          dvdrho_out=dvdrho_1d,iwarn_dvdrho=iwarn_dv)

     if(ierr.ne.0) then
        write(lunerr,*) ' ?eq_flxint_arr2d: xplasma_rho_zonint call failed.'

     else

        if(iwant.eq.3) then 
           result(1,1:inrhoi-iadj)= &
                result(1,1:inrhoi-iadj)*dvol_1d(1:inrhoi-iadj)

        else if(iwant.eq.5) then
           result(1,1:inrhoi)=result(1,1:inrhoi)*dvdrho_1d(1:inrhoi)

        endif
     endif
     deallocate(dvol_1d,dvdrho_1d)

  else
     !  integrate over 2d set of zones or segments of surfaces

     inthzons = inthi-1

     if(inumrho.lt.inrhoi-iadj) then
        write(lunerr,*) ' ?eq_flxint_arr2d: "inumrho" dimension too small.'
        write(lunerr,*) '  found: ',inumrho,' need at least: ',inrhoi-iadj
        ierr=1
     endif

     if(inumchi.lt.inthzons) then
        write(lunerr,*) ' ?eq_flxint_arr2d: "inumchi" dimension too small.'
        write(lunerr,*) '  found: ',inumchi,' need at least: ',inthzons
        ierr=1
     endif
     if(ierr.ne.0) return

     allocate(dvol_2d(inthzons,inrhoi-1),dvdrho_2d(inthzons,inrhoi))

     call xplasma_2d_zonint(s,idi,zname,result(1:inthzons,1:inrhoi-iadj),ierr,&
          dvol_out=dvol_2d,iwarn_dvol=iwarn_vol, &
          dvdrho_out=dvdrho_2d,iwarn_dvdrho=iwarn_dv)

     if(ierr.ne.0) then
        write(lunerr,*) ' ?eq_flxint_arr2d: xplasma_2d_zonint call failed.'

     else

        if(iwant.eq.3) then 
           result(1:inthzons,1:inrhoi-iadj)= &
                result(1:inthzons,1:inrhoi-iadj)*dvol_2d(1:inthzons,1:inrhoi-iadj)

        else if(iwant.eq.5) then
           result(1:inthzons,1:inrhoi)= &
                result(1:inthzons,1:inrhoi)*dvdrho_2d(1:inthzons,1:inrhoi)

        endif
     endif
     deallocate(dvol_2d,dvdrho_2d)

  endif

  !  (3)  cleanup

  call xoi_author_set(iertmp)
  call xplasma_remove_item(s,id_tmp,iertmp)
  call xoi_author_clear(iertmp)

end subroutine eq_flxint_arr2d
