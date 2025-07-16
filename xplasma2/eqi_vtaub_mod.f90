module eqi_vtaub_mod

  !  module to set up splines for faster tau_bounce integral evaluation
  !  adaptive integrator needed due to sqrt singularity at bounce point

  implicit NONE
  SAVE

  integer, parameter :: imax=1024

  real*8, parameter :: bcxmin=0.0d0, bcxmax=0.0d0
  integer, parameter :: ibcxmin=-1, ibcxmax=-1

  !  support independent spline sets for multiple threads (up to imax...)

  type :: vtaub_data
     !  assume f95 OK -- allocatable elements

     integer :: nth = 0
     real*8, dimension(:), allocatable :: th,wk
     real*8, dimension(:,:), allocatable :: fmodb,fbpol,fdldth
  end type vtaub_data

  type (vtaub_data) :: vta(imax)

  real*8, parameter :: ZPI = 3.1415926535897931D+00

  contains

    subroutine elem_alloc(indx,inth,ierr)

      !  check that index is in range 0 to imax; allocate spline arrays

      integer, intent(in) :: indx    ! index of type element to use...
      integer, intent(in) :: inth    ! no. of theta points -pi:pi inclusive
      integer, intent(out) :: ierr   ! 0=OK

      !----------------
      integer :: ith,inthi
      !----------------

      inthi=max(41,min(1001,inth))  ! sanity check

      if((indx.lt.1).or.(indx.gt.imax)) then
         ierr=1
      else
         ierr = 0
         if(allocated(vta(indx)%th)) then
            deallocate(vta(indx)%th)
            deallocate(vta(indx)%wk)
            deallocate(vta(indx)%fmodb)
            deallocate(vta(indx)%fbpol)
            deallocate(vta(indx)%fdldth)
         endif
         vta(indx)%nth = inthi
         allocate(vta(indx)%th(inthi))
         allocate(vta(indx)%wk(inthi))
         allocate(vta(indx)%fmodb(4,inthi))
         allocate(vta(indx)%fbpol(4,inthi))
         allocate(vta(indx)%fdldth(4,inthi))
         do ith=1,inthi
            vta(indx)%th(ith)= -ZPI + (ith-1)*2*ZPI/(inthi-1)
         enddo
      endif

    end subroutine elem_alloc

    subroutine elem_init(indx,fmodb,fbpol,fdldth,ierr)

      !  fill spline arrays for indicated element

      integer, intent(in) :: indx           ! element index
      real*8, dimension(:), intent(in) :: fmodb  ! |B| vs. theta
      real*8, dimension(:), intent(in) :: fbpol  ! |Bpol| vs. theta
      real*8, dimension(:), intent(in) :: fdldth ! |dl/dth| vs. theta
      integer, intent(out) :: ierr          ! exit status, 0=OK

      !-------------------
      integer :: isize_min,isize_max,inth,idum,iertmp
      !-------------------

      ierr=0
      if((indx.lt.1).or.(indx.gt.imax)) then
         ierr=1
      else if(.not.allocated(vta(indx)%th)) then
         ierr=2
      else
         isize_min=size(fmodb)
         isize_max=isize_min

         isize_min=min(isize_min,size(fbpol))
         isize_max=max(isize_max,size(fbpol))

         isize_min=min(isize_min,size(fdldth))
         isize_max=max(isize_max,size(fdldth))

         if((isize_min.ne.vta(indx)%nth).or.(isize_max.ne.vta(indx)%nth)) then
            ierr=3
         endif
      endif

      if(ierr.ne.0) return

      inth=vta(indx)%nth

      vta(indx)%fmodb(1,1:inth) = fmodb
      vta(indx)%fbpol(1,1:inth) = fbpol
      vta(indx)%fdldth(1,1:inth) = fdldth

      call r8cspline(vta(indx)%th,inth, &
           vta(indx)%fmodb, &
           ibcxmin,bcxmin,ibcxmax,bcxmax, &
           vta(indx)%wk,inth,idum,iertmp)
      if(iertmp.ne.0) ierr=4

      call r8cspline(vta(indx)%th,inth, &
           vta(indx)%fbpol, &
           ibcxmin,bcxmin,ibcxmax,bcxmax, &
           vta(indx)%wk,inth,idum,iertmp)
      if(iertmp.ne.0) ierr=4

      call r8cspline(vta(indx)%th,inth, &
           vta(indx)%fdldth, &
           ibcxmin,bcxmin,ibcxmax,bcxmax, &
           vta(indx)%wk,inth,idum,iertmp)
      if(iertmp.ne.0) ierr=4
      
    end subroutine elem_init
      
    subroutine errstr(ierr,str)

      !  return a text string describing an error

      integer, intent(in) :: ierr
      character*(*), intent(out) :: str

      if(ierr.eq.1) then
         str='eqi_vtaub_mod: element index not between 1 and imax.'

      else if(ierr.eq.2) then
         str='eqi_vtaub_mod: element not allocated.'

      else if(ierr.eq.3) then
         str='eqi_vtaub_mod: size of passed data does not match element grid size.'

      else if(ierr.eq.4) then
         str='eqi_vtaub_mod: spline setup error.'

      else if(ierr.ne.0) then
         str='eqi_vtaub_mod: error (unspecified).'

      else
         str=' '
      endif

    end subroutine errstr  

end module eqi_vtaub_mod
