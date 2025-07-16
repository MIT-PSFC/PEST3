subroutine geqdsk_index1_ntimes(ilun,filename,ntimes)

  use geqdsk_index_mod
  implicit NONE

  integer, intent(in) :: ilun  ! fortran LUN for file input
  character*(*), intent(in) :: filename  ! file to read (index of G-eqdsks...)
  integer, intent(out) :: ntimes  ! no. of times indexed in the file

  !  get number of times; return ntimes=0 if there is an error...

  !---------------
  integer :: ilun_msg,ierr
  !---------------

  if(.not.geqdsk_index1_init) then
     call geqdsk_index_preset(geqdsk_index1)
     geqdsk_index1_init=.TRUE.
  endif

  call geqdsk_index_init(ilun,filename,geqdsk_index1,ierr)
  if(ierr.ne.0) then
     call geqdsk_lunmsg_get(ilun_msg)
     write(ilun_msg,*) ' ?geqdsk_index1_ntimes: error reading file: ',trim(filename)
     call geqdsk_index_msg(geqdsk_index1,ierr,ilun_msg)
     ntimes=0
  else
     ntimes=geqdsk_index1%ntimes
  endif

end subroutine geqdsk_index1_ntimes

subroutine geqdsk_index1_path(ilun,filename,path)

  use geqdsk_index_mod
  implicit NONE

  integer, intent(in) :: ilun  ! fortran LUN for file input
  character*(*), intent(in) :: filename  ! file to read (index of G-eqdsks...)
  character*(*), intent(out) :: path     ! path to files (blank if $cwd)

  !  get number of times; return path=0 if there is an error...
  !---------------
  integer :: ilun_msg,ierr
  !---------------

  if(.not.geqdsk_index1_init) then
     call geqdsk_index_preset(geqdsk_index1)
     geqdsk_index1_init=.TRUE.
  endif

  call geqdsk_index_init(ilun,filename,geqdsk_index1,ierr)
  if(ierr.ne.0) then
     call geqdsk_lunmsg_get(ilun_msg)
     write(ilun_msg,*) ' ?geqdsk_index1_path: error reading file: '&
          &,trim(filename)
     call geqdsk_index_msg(geqdsk_index1,ierr,ilun_msg)
     path='?'
  else
     path=geqdsk_index1%path
  endif

end subroutine geqdsk_index1_path

subroutine geqdsk_index1_gtime_r4(ilun,filename,itime,gtime,ierr)
  implicit NONE

  integer, intent(in) :: ilun  ! fortran LUN for file input
  character*(*), intent(in) :: filename  ! file to read (index of G
  !-eqdsks...)
  integer, intent(in) :: itime           ! time index
  real, intent(out) :: gtime             ! time value returned (r4)
  integer, intent(out) :: ierr           ! completion code (0=OK)

  real*8 :: gtime_r8

  call geqdsk_index1_gtime(ilun,filename,itime,gtime_r8,ierr)
  gtime=gtime_r8

end subroutine geqdsk_index1_gtime_r4

subroutine geqdsk_index1_gtime(ilun,filename,itime,gtime,ierr)

  ! return the itime'th time in the index file (filename)

  use geqdsk_index_mod
  implicit NONE

  integer, intent(in) :: ilun  ! fortran LUN for file input
  character*(*), intent(in) :: filename  ! file to read (index of G
  !-eqdsks...)
  integer, intent(in) :: itime           ! time index
  real*8, intent(out) :: gtime           ! time value returned
  integer, intent(out) :: ierr           ! completion code (0=OK)

  !  get time associated with index "itime" in the index file "filename".

  !---------------
  integer :: ilun_msg
  !---------------

  if(.not.geqdsk_index1_init) then
     call geqdsk_index_preset(geqdsk_index1)
     geqdsk_index1_init=.TRUE.
  endif

  call geqdsk_index_init(ilun,filename,geqdsk_index1,ierr)

  gtime=0
  if(ierr.ne.0) then
     call geqdsk_lunmsg_get(ilun_msg)
     write(ilun_msg,*) ' ?geqdsk_index1_gtime: error reading file: ',trim(filename)
     call geqdsk_index_msg(geqdsk_index1,ierr,ilun_msg)
  else if((itime.le.0).or.(itime.gt.geqdsk_index1%ntimes)) then
     call geqdsk_lunmsg_get(ilun_msg)
     write(ilun_msg,*) ' ?geqdsk_index1_gtime: ',trim(filename)
     write(ilun_msg,*) '  time index: ',itime,' out of file index range: ', &
          '1 to ',geqdsk_index1%ntimes
  else
     gtime=geqdsk_index1%times(itime)
  endif

end subroutine geqdsk_index1_gtime

subroutine geqdsk_index1_gfile(ilun,filename,itime,gfile,ierr)

  ! return the itime'th gfile filename in the index file (filename)

  use geqdsk_index_mod
  implicit NONE

  integer, intent(in) :: ilun  ! fortran LUN for file input
  character*(*), intent(in) :: filename  ! file to read (index of G-eqdsks...)
  integer, intent(in) :: itime           ! time index
  character*(*), intent(out) :: gfile    ! G-eqdsk filename returned
  integer, intent(out) :: ierr           ! completion code (0=OK)

  !  get number of times; return ntimes=0 if there is an error...

  !---------------
  integer :: ilun_msg
  !---------------

  if(.not.geqdsk_index1_init) then
     call geqdsk_index_preset(geqdsk_index1)
     geqdsk_index1_init=.TRUE.
  endif

  call geqdsk_index_init(ilun,filename,geqdsk_index1,ierr)

  gfile=' '
  if(ierr.ne.0) then
     call geqdsk_lunmsg_get(ilun_msg)
     write(ilun_msg,*) ' ?geqdsk_index1_gfile: error reading file: ',trim(filename)
     call geqdsk_index_msg(geqdsk_index1,ierr,ilun_msg)
  else if((itime.le.0).or.(itime.gt.geqdsk_index1%ntimes)) then
     call geqdsk_lunmsg_get(ilun_msg)
     write(ilun_msg,*) ' ?geqdsk_index1_gfile: ',trim(filename)
     write(ilun_msg,*) '  time index: ',itime,' out of file index range: ', &
          '1 to ',geqdsk_index1%ntimes
  else
     gfile=geqdsk_index1%efiles(itime)
  endif

end subroutine geqdsk_index1_gfile

subroutine geqdsk_findex(fstring,iend,indx)

  ! given a string of the form "EFIT_INDEX:<file-path>(<integer-index>)
  ! return the start and stop string indices of <file-path>, and decode
  ! the integer index inside the parentheses; return iend=indx=0 if an
  ! error occurs; write no messages.

  implicit NONE
  character*(*), intent(in) :: fstring  ! input string
  integer,intent(out) :: iend           ! output loc. end of <file-path> string
  integer,intent(out) :: indx           ! output decoded <integer-index> value

  character*10 zstr
  integer i,j,ilparen,irparen,ili,istat

  iend=0
  indx=0

  zstr=fstring(1:10)
  call uupper(zstr)
  if(zstr//fstring(11:11).ne.'EFIT_INDEX:') then
     return   ! no lead string
  endif

  ilparen=index(fstring,'(')
  irparen=index(fstring,')')
  if(min(ilparen,irparen).le.0) then
     return   ! no index value or missing parentheses
  endif

  j=11
  zstr=' '
  do i=irparen-1,ilparen+1,-1
     if(fstring(i:i).ne.' ') then
        j=j-1
        if(j.eq.0) exit
        zstr(j:j)=fstring(i:i)
     endif
  enddo
  if((j.eq.11).or.(j.eq.0)) then
     return   ! empty or too long integer string
  endif

  read(zstr,'(i10)',iostat=istat) indx
  if(istat.ne.0) then
     indx=0
     return   ! decode failure
  endif

  ! if we made it this far, the index decode was successful
  ! although it is still possible for the value to be out of range
  ! with respect to the contents of the named file...

  ! find end of filename and exit...

  do i=ilparen-1,11,-1
     if(fstring(i:i).ne.' ') then
        iend=i
        exit
     endif
  enddo

end subroutine geqdsk_findex
