!  hand maintained R4 interfaces ... for eqm_rhofun & eqm_frhochi
!                                ... and eq_hrhochi ... etc ...
c----------------------------------------------
      subroutine r4_eqm_rhofun(iorder,id_axis,zlbl,r4_zdata,
     >   ibc1,r4_zbc1,ibc2,r4_zbc2,
     >   id_rhofun,ierr)
c
c  REAL (usually real*4) interface
cc  establish a function of "rho" -- axis to plasma bdy
c
c  input:
      IMPLICIT NONE
c
      integer iorder                    ! interpolation order
c
c     for all of the following, (nrho) values are to be provided.
c
c  iorder = 0 -- piecewise linear
c  iorder = 1 -- use Hermite piecewise cubics with Akima method for
c     setting the derivative at the grid points (modified at the bdys)
c  iorder = 2 or 3 -- use a cubic spline
c
c  iorder = 1 -- once differentiable
c  iorder = 2 or 3 -- twice differentiable (spline)
c
      integer id_axis                   ! axis id -- must be rho or akin to rho
c
c  if id_axis is not rho, then the size of the object is not nrho but the
c  size of the given axis.  The given axis must be "kin" to rho, i.e. it
c  must map the same space, but perhaps with a larger or smaller number of
c  points differently distributed.
c
      character*(*) zlbl                ! function name (label)
      REAL r4_zdata(*)                  ! data (at least `nrho' words)
c
      integer ibc1                      ! bc @ rho_axis
      real r4_zbc1                      ! bc parameter
      integer ibc2                      ! bc @ rho_bdy
      real r4_zbc2                      ! bc parameter
c
c                =0:  "standard" boundary condition:  "not a knot" for
c                     splines, Akima for piecewise Hermite.
c                =1:  assigned df/drho bc, --> df/drho=zbc1|2
c                =2:  assigned d2f/drho2 bc (splines only), d2f/drho2=zbc1|2
c
c  for iorder.le.0, bc controls are ignored.
c
c  output:
c
      integer id_rhofun                 ! function id, returned
      integer ierr                      ! completion code, 0=OK
c
c-----------------------------
c
      real*8, dimension(:), allocatable :: zdata
      real*8 zbc1,zbc2
c
      integer ixsize
c
c-----------------------------
c
      call eq_ngrid(id_axis,ixsize)
c
      if(ixsize.gt.0) then
         allocate(zdata(ixsize))
         zdata=r4_zdata(1:ixsize)
      else
         allocate(zdata(1))
         zdata=0  ! there will be an error...
      endif
c
      zbc1=r4_zbc1
      zbc2=r4_zbc2
c
      call eqm_rhofun(iorder,id_axis,zlbl,zdata,
     >   ibc1,zbc1,ibc2,zbc2,
     >   id_rhofun,ierr)
c
      deallocate(zdata)
c
      return
      end
c----------------------------------------------
      subroutine r4_eqm_f1d(iorder,id_axis,zlbl,r4_zdata,
     >   ibc1,r4_zbc1,ibc2,r4_zbc2,
     >   id_f1d,ierr)
c
c  REAL (usually real*4) interface
cc  establish a function of "rho" -- axis to plasma bdy
c
c  input:
      IMPLICIT NONE
c
      integer iorder                    ! interpolation order
c
c     for all of the following, (nrho) values are to be provided.
c
c  iorder = 0 -- piecewise linear
c  iorder = 1 -- use Hermite piecewise cubics with Akima method for
c     setting the derivative at the grid points (modified at the bdys)
c  iorder = 2 or 3 -- use a cubic spline
c
c  iorder = 1 -- once differentiable
c  iorder = 2 or 3 -- twice differentiable (spline)
c
      integer id_axis                   ! axis id -- must be rho or akin to rho
c
c  if id_axis is not rho, then the size of the object is not nrho but the
c  size of the given axis.  The given axis must be "kin" to rho, i.e. it
c  must map the same space, but perhaps with a larger or smaller number of
c  points differently distributed.
c
      character*(*) zlbl                ! function name (label)
      REAL r4_zdata(*)                  ! data (at least `nrho' words)
c
      integer ibc1                      ! bc @ rho_axis
      real r4_zbc1                      ! bc parameter
      integer ibc2                      ! bc @ rho_bdy
      real r4_zbc2                      ! bc parameter
c
c                =0:  "standard" boundary condition:  "not a knot" for
c                     splines, Akima for piecewise Hermite.
c                =1:  assigned df/drho bc, --> df/drho=zbc1|2
c                =2:  assigned d2f/drho2 bc (splines only), d2f/drho2=zbc1|2
c
c  for iorder.le.0, bc controls are ignored.
c
c  output:
c
      integer id_f1d                 ! function id, returned
      integer ierr                      ! completion code, 0=OK
c
c-----------------------------
c
      real*8, dimension(:), allocatable :: zdata
      real*8 zbc1,zbc2
c
      integer ixsize
c
c-----------------------------
c
      call eq_ngrid(id_axis,ixsize)
c
      if(ixsize.gt.0) then
         allocate(zdata(ixsize))
         zdata=r4_zdata(1:ixsize)
      else
         allocate(zdata(1))
         zdata=0  ! there will be an error...
      endif
c
      zbc1=r4_zbc1
      zbc2=r4_zbc2
c
      call eqm_f1d(iorder,id_axis,zlbl,zdata,
     >   ibc1,zbc1,ibc2,zbc2,
     >   id_f1d,ierr)
c
      deallocate(zdata)
c
      return
      end
c----------------------------------------------------------
      subroutine r4_eqm_frhochi(iorder,id_ax1,id_ax2,zlbl,r4_zdata,id1,
     >   ibcrho0,r4_zbcrho0,ibcrho1,r4_zbcrho1,id_fun,ierr)
c
c  REAL (usually real*4) interface to eqm_frhochi
c
c  establish a function of "rho & chi" -- axis to plasma bdy
c
      use xplasma_obj_instance
c
      IMPLICIT NONE
c
c  input:
      integer iorder                    ! desired interpolation order
c
c  [in these comments nrho & nchi refer to the axes of the data passed]
c  iorder=-1: (nrho-1)x(nchi-1) pts supplied -- binned data, with
c        f(i,j) constant in [rho(i),rho(i+1)]x[chi(j),chi(j+1)]
c
c  iorder=0:  bilinear interpolation on data on (nrho)x(nchi) grid
c  iorder=1:  Akima Hermite interpolation on data on(nrho)x(nchi) grid
c  iorder=2 or 3:  Spline interpolation on data on (nrho)x(nchi) grid
c
      integer id_ax1                    ! axis id -- 1st dimension of zdata
      integer id_ax2                    ! axis id -- 2nd dimension of zdata
c
c  id_ax1 or id_ax2 must give a rho coordinate; the other coordinate must
c  be a chi (poloidal angle) coordinate.  Regardless of the ordering of
c  zdata, the module is data is stored with chi = 1st dimension, rho =
c  2nd dimension.
c
      character*(*) zlbl                ! function name (label)
      integer id1                       ! size of 1st dimension of zdata
c
      REAL r4_zdata(id1,*)              ! data (at least nrho*nchi words)
c
c  id1 must be .ge. size of axis (id_ax1)
c
c  bdy conditions are enforced only for iorder.ge.1
c
      integer ibcrho0,ibcrho1           ! BC's at extrema
c
c  BC controls are as in "pspline":  0 -- default, 1 -- fix df/drho
c                                    2 -- fix d2f/drho2
c
      REAL r4_zbcrho0(*)                ! if ibc>0:  df/drho @ rho(1)
      REAL r4_zbcrho1(*)                ! ditto:  df/drho @ rho(nrho)
c
c  output:
c
      integer id_fun                    ! function id code
      integer ierr                      ! completion code, 0=OK
c
c-----------------------------
c  local copy
c
      real*8, dimension(:,:), allocatable :: zdata
      real*8, dimension(:), allocatable :: zbcrho0,zbcrho1
c
      integer i,in1,in2,inrho,inchi,iflag,icoord
c
c-----------------------------
c
c  find dimensions...
c
      call eq_ngrid(id_ax1,in1)
      call eq_ngrid(id_ax2,in2)

      iflag=0
      if((in1.eq.0).or.(in2.eq.0)) iflag=1

      if(iflag.eq.0) then
c
c  axes look OK; which one is "rho"?
c
         iflag=1
c
         call xplasma_grid_info(s,id_ax1,ierr, coord=icoord)
         if(icoord.eq.xplasma_rho_coord) then
            inrho=in1
            inchi=in2
            iflag=0

         else
            call xplasma_grid_info(s,id_ax2,ierr, coord=icoord)
            if(icoord.eq.xplasma_rho_coord) then
               inrho=in2
               inchi=in1
               iflag=0
            endif
         endif
      endif

      if(iflag.eq.1) then
         allocate(zbcrho0(1),zbcrho1(1),zdata(id1,1))
         zbcrho0=0; zbcrho1=0; zdata=0
      else
         allocate(zbcrho0(inchi))
         allocate(zbcrho1(inchi))
         zbcrho0=0
         zbcrho1=0
         if(ibcrho0.gt.0) zbcrho0=r4_zbcrho0(1:inchi)
         if(ibcrho1.gt.0) zbcrho1=r4_zbcrho1(1:inchi)
c
         allocate(zdata(id1,in2))
         zdata=0
         do i=1,in2
            zdata(1:in1,i)=r4_zdata(1:in1,i)
         enddo
c
      endif
c
      call eqm_frhochi(iorder,id_ax1,id_ax2,zlbl,zdata,id1,
     >     ibcrho0,zbcrho0,ibcrho1,zbcrho1,id_fun,ierr)
c
      deallocate(zdata)
      deallocate(zbcrho1)
      deallocate(zbcrho0)
c
      return
      end
c------------------------------------------------------
      subroutine r4_eq_hrhochi(ivec,r4_zrho,r4_zchi,nlist,iflist,
     >   ictwant,ivecd,r4_zans,ierr)
c
c  hand maintained real*4 interface to eq_hrhochi
c  (dmc 25 Aug 2000)
c
      implicit NONE
c
c  input:
c
      integer ivec                      ! input vector dimension
      !  ivec.lt.0 to indicate CW-reversed chi

      REAL r4_zrho(abs(ivec))
      REAL r4_zchi(abs(ivec))
      integer nlist                     ! number of functions
      integer iflist(nlist)             ! functions to evaluate
c
      integer ictwant(6)                ! output selector (pspline style)
c
c  each element of ictwant(1:6) should be 0 or 1
c  use 1 to indicate a desired item:
c
c   ictwant(1) -- function values            ...1 if wanted, 0 otherwise
c   ictwant(2) -- functions' df/drho values
c   ictwant(3) -- functions' df/dchi values
c   ictwant(4) -- functions' d2f/drho2 values
c   ictwant(5) -- functions' d2f/dchi2 values
c   ictwant(6) -- functions' d2f/drhodchi values
c
      integer ivecd                     ! output vector dimension
c
c  output:
c
      REAL r4_zans(ivecd,*)              ! evaluation results
      integer ierr                      ! completion code (0=OK)
c
c  ** let ictnum = # of non-zero elements of ictwant.
c         (ictnum=sum(ictwant)).  Then...
c  zans(1:ivec,1)   = value or derivative for function id iflist(1)
c                     indicated by first non-zero value of ictwant(...)
c  zans(1:ivec,2)   = value or derivative for function id iflist(1)
c                     indicated by second non-zero value of ictwant(...)
c    ...
c  zans(1:ivec,ictnum)   = value or derivative for function id iflist(1)
c                     indicated by last non-zero value of ictwant(...)
c
c  or generally, for 1 .le. k .le. ictnum,
c    ...
c  zans(1:ivec,(j-1)*ictnum+k) = value or derivative for function
c                                id iflist(j), corresponding to k'th
c                                non-zero element of ictwant(...)
c
c--------------------------------------------------
c
      REAL*8 zrho(abs(ivec))                 ! argument -- where to evaluate
      REAL*8 zchi(abs(ivec))                 ! argument -- where to evaluate
c
      REAL*8, dimension(:,:), allocatable :: zans
c
      integer ictnum,i
c
c----------------------------------
c
c  copy to REAL*8 argument arrays
c
      zrho=r4_zrho
      zchi=r4_zchi
c
c  allocate REAL*8 output buffer
c
      ictnum=sum(ictwant)
      allocate(zans(ivecd,max(1,ictnum)*nlist))
c
c  call eq_hrhochi using REAL*8 arguments
c
      call eq_hrhochi(ivec,zrho,zchi,nlist,iflist,ictwant,
     >   ivecd,zans,ierr)
c
c  copy answer back
c
      if(ierr.eq.0) then
c
         do i=1,max(1,ictnum)*nlist
            r4_zans(1:abs(ivec),i)=zans(1:abs(ivec),i)
         enddo
c
      endif
c
c  deallocate buffer
c
      deallocate(zans)
c
      return
      end
c------------------------------------------------------
      subroutine r4_eq_hRZ(ivec,r4_zR,r4_zZ,nlist,iflist,
     >   ictwant,ivecd,r4_zans,ierr)
c
c  hand maintained real*4 interface to eq_hRZ
c  (dmc 25 Aug 2000)
c
      implicit NONE
c
c  input:
c
      integer ivec                      ! input vector dimension
      REAL r4_zR(ivec)
      REAL r4_zZ(ivec)
      integer nlist                     ! number of functions
      integer iflist(nlist)             ! functions to evaluate
c
      integer ictwant(6)                ! output selector (pspline style)
c
c  each element of ictwant(1:6) should be 0 or 1
c  use 1 to indicate a desired item:
c
c   ictwant(1) -- function values            ...1 if wanted, 0 otherwise
c   ictwant(2) -- functions' df/dR values
c   ictwant(3) -- functions' df/dZ values
c   ictwant(4) -- functions' d2f/dR2 values
c   ictwant(5) -- functions' d2f/dZ2 values
c   ictwant(6) -- functions' d2f/dRdZ values
c
      integer ivecd                     ! output vector dimension
c
c  output:
c
      REAL r4_zans(ivecd,*)             ! evaluation results
      integer ierr                      ! completion code (0=OK)
c
c  ** let ictnum = # of non-zero elements of ictwant.
c         (ictnum=sum(ictwant)).  Then...
c  zans(1:ivec,1)   = value or derivative for function id iflist(1)
c                     indicated by first non-zero value of ictwant(...)
c  zans(1:ivec,2)   = value or derivative for function id iflist(1)
c                     indicated by second non-zero value of ictwant(...)
c    ...
c  zans(1:ivec,ictnum)   = value or derivative for function id iflist(1)
c                     indicated by last non-zero value of ictwant(...)
c
c  or generally, for 1 .le. k .le. ictnum,
c    ...
c  zans(1:ivec,(j-1)*ictnum+k) = value or derivative for function
c                                id iflist(j), corresponding to k'th
c                                non-zero element of ictwant(...)
c
c--------------------------------------------------
c
      REAL*8 zR(ivec)                   ! argument -- where to evaluate
      REAL*8 zZ(ivec)                   ! argument -- where to evaluate
c
      REAL*8, dimension(:,:), allocatable :: zans
c
      integer ictnum,i
c
c----------------------------------
c
c  copy to REAL*8 argument arrays
c
      zR=r4_zR
      zZ=r4_zZ
c
c  allocate REAL*8 output buffer
c
      ictnum=sum(ictwant)
      allocate(zans(ivecd,max(1,ictnum)*nlist))
c
c  call eq_hRZ using REAL*8 arguments
c
      call eq_hRZ(ivec,zR,zZ,nlist,iflist,ictwant,
     >   ivecd,zans,ierr)
c
c  copy answer back
c
      if(ierr.eq.0) then
c
         do i=1,max(1,ictnum)*nlist
            r4_zans(1:ivec,i)=zans(1:ivec,i)
         enddo
c
      endif
c
c  deallocate buffer
c
      deallocate(zans)
c
      return
      end
