subroutine eqm_rzmag(Rarr,Zarr,id1,id2,idrho,id_R,id_Z,ierr)
 
  use xplasma_obj_instance
  use eq_module
  implicit NONE

  !  establish R(chi,rho) and Z(chi,rho) bicubic splines.
  !  periodic in chi, not-a-knot BC @ rho bdys
  !    (for more control of rho BCs use eqm_rzmagb)

  !  input arguments:

  integer id1,id2            ! R,Z array dimensions
  real*8 Rarr(id1,id2)       ! R array
  real*8 Zarr(id1,id2)       ! Z array

  integer idrho              ! =1:  id1 is "rho" dimension; =2: id2 is.
                             ! =-1: id1 is "rho"; use finite diff BC
                             ! =-2: id2 is "rho"; use finite diff BC

  !  output arguments:

  integer id_R               ! xplasma id for "R" bicubic spline (returned)
  integer id_Z               ! xplasma id for "Z" bicubic spline (returned)

  integer ierr               ! completion code (0 = OK)
  
  !----------------------------------------------------------------

  real*8, dimension(:), allocatable :: zrbc0,zrbc1,zzbc0,zzbc1,zrho
  real*8 :: zdrho0,zdrho1
  integer ibc  ! will use 1st derivative BC
  integer i,inrho,inchi,id_rho,id_chi,inchi_err

  !----------------------------------------------------------------

  call eqm_rzmag_check(id1,id2,abs(idrho),id_rho,id_chi,inrho,inchi,ierr)
  if(ierr.ne.0) return
  
  allocate(zrho(inrho))
  call xplasma_grid(s,id_rho,zrho,ierr)
  if(ierr.ne.0) then
     write(lunerr,*) ' ?eqm_rzmag: "__RHO" grid fetch failed.'
     call xplasma_error(s,ierr,lunerr)
     return
  endif

  ! use finite difference formula for 1st derivative BCs.  This seems to
  ! be safest for preserving the 2x2 Jacobian sign...

  zdrho0=zrho(2)-zrho(1)
  zdrho1=zrho(inrho)-zrho(inrho-1)

  if(idrho.eq.-1) then
     allocate(zrbc0(id2),zrbc1(id2),zzbc0(id2),zzbc1(id2))
     if(id2.lt.inchi) then
        ierr=1
        inchi_err=id2
     else
        do i=1,id2
           if(i.le.inchi) then
              zrbc0(i)=(Rarr(2,i)-Rarr(1,i))/zdrho0
              zrbc1(i)=(Rarr(inrho,i)-Rarr(inrho-1,i))/zdrho1
              zzbc0(i)=(Zarr(2,i)-Zarr(1,i))/zdrho0
              zzbc1(i)=(Zarr(inrho,i)-Zarr(inrho-1,i))/zdrho1
           else
              zrbc0(i)=0
              zrbc1(i)=0
              zzbc0(i)=0
              zzbc1(i)=0
           endif
        enddo
     endif

  else if(idrho.eq.-2) then
     allocate(zrbc0(id1),zrbc1(id1),zzbc0(id1),zzbc1(id1))
     if(id1.lt.inchi) then
        ierr=1
        inchi_err=id1
     else
        do i=1,id1
           if(i.le.inchi) then
              zrbc0(i)=(Rarr(i,2)-Rarr(i,1))/zdrho0
              zrbc1(i)=(Rarr(i,inrho)-Rarr(i,inrho-1))/zdrho1
              zzbc0(i)=(Zarr(i,2)-Zarr(i,1))/zdrho0
              zzbc1(i)=(Zarr(i,inrho)-Zarr(i,inrho-1))/zdrho1
           else
              zrbc0(i)=0
              zrbc1(i)=0
              zzbc0(i)=0
              zzbc1(i)=0
           endif
        enddo
     endif

  else
     allocate(zrbc0(inchi),zrbc1(inchi),zzbc0(inchi),zzbc1(inchi))
     zrbc0=0
     zrbc1=0
     zzbc0=0
     zzbc1=0
  endif

  if(ierr.eq.1) then
     write(lunerr,*) ' ?eqm_rzmag: too few chi points in passed array.'
     write(lunerr,*) '  needed: ',inchi,' got: ',inchi_err
     write(lunerr,*) '  dimension error in poloidal angle coordinate.'

  else

     if(abs(idrho).lt.0) then
        ibc=1
     else
        ibc=0
     endif

     call eqm_rzmagb(Rarr,Zarr,id1,id2,idrho,  &
          ibc,zrbc0,ibc,zrbc1,  &
          ibc,zzbc0,ibc,zzbc1,  &
          id_R,id_Z,ierr)

  endif

  deallocate(zrbc0,zrbc1,zzbc0,zzbc1)

end subroutine eqm_rzmag
