program proclist

  ! write file containing list of processors and processor ID #s in current job
  ! filename = 1st command line argument, which must be nonblank
  !    if there are more than one command line arguments, these are ignored.

  implicit NONE

#ifdef __MPI
  include 'mpif.h'
#endif

  ! local communicator

#ifdef __MPI
  integer :: nbi_comm = MPI_COMM_WORLD
#endif

  !-------------------------------------
  character*180 filenam,hostnam,hostnam1,proc_name
  integer :: nargs,istat,myid,numprocs,ierr
  integer :: ii,jj,mywid,maxwid,iresult,idum
  integer :: lunout = 77
  integer, dimension(:), allocatable :: iwrite

  integer, dimension(:,:), allocatable :: intbuf,intbufr
  !-------------------------------------
  ! executable code

#ifdef __MPI

  call MPI_INIT(ierr)
  if(ierr.ne.0) call errmsg_exit(' MPI_INIT failure! ')

  call MPI_COMM_RANK(nbi_comm,myid,ierr);   call ckerr('MPI_COMM_RANK')
  call MPI_COMM_SIZE(nbi_comm,numprocs,ierr);   call ckerr('MPI_COMM_SIZE')
  if(myid.eq.0) then
     write(0,*) ' mpi_proclist:  numprocs = ',numprocs
  endif

#else

  ierr=0

  write(0,*) ' proclist: serial mode.'
  myid=0
  numprocs=0

#endif

  if (allocated(iwrite))deallocate(iwrite)
  allocate(iwrite(0:numprocs-1))
  iwrite(:)=0

  if(myid.eq.0) then

     call get_arg_count(nargs)
     if(nargs.le.0) then

        write(0,*) ' ?filename argument missing.'
        write(0,*) ' ?syntax:  [mpi_]proclist  <filename>'
        call exit(1)

     endif

     call get_arg(1,filenam)
     open(unit=lunout,file=filenam,status='unknown',iostat=istat)

     if(istat.ne.0) then

        write(0,*) ' ?filename open failure: "',trim(filenam),'".'
        write(0,*) ' ?syntax:  [mpi_]proclist  <filename>'
        call exit(1)

     endif
  endif

#ifdef __MPI

  call MPI_GET_PROCESSOR_NAME(proc_name,iresult,ierr)
  call ckerr('MPI_GET_PROCESSOR_NAME')

  mywid = len(trim(proc_name))

  call MPI_ALLREDUCE(mywid,maxwid,1, &
       MPI_INTEGER,MPI_MAX,nbi_comm,ierr)
  call ckerr('MPI_ALLREDUCE(mywid)')

  allocate(intbuf(maxwid,0:numprocs-1)); intbuf = 0
  allocate(intbufr(maxwid,0:numprocs-1))

  do ii=1,mywid
     intbuf(ii,myid) = ichar(proc_name(ii:ii))
  enddo

  call MPI_REDUCE(intbuf,intbufr, numprocs*maxwid, &
       MPI_INTEGER,MPI_SUM,0,nbi_comm,ierr)
  call ckerr('MPI_REDUCE(intbuf)')

  if(myid.eq.0) then
     hostnam1=' '
     do jj=0,numprocs-1
        hostnam=' '
        do ii=1,maxwid
           if(intbufr(ii,jj).eq.0) exit
           hostnam(ii:ii) = char(intbufr(ii,jj))
        enddo
!       if (trim(hostnam).ne.trim(hostnam1))iwrite(jj)=1
        if (trim(hostnam).ne.trim(hostnam1))iwrite(jj)=jj
        if (jj.eq.0) iwrite(jj)=0
        hostnam1=hostnam
        write(lunout,'(1x,i8,1x,a,1x,i8)') jj, trim(hostnam), iwrite(jj)
     enddo
  endif

#else
  call sget_host(hostnam)

  idum = 0
  write(lunout,'(1x,i8,1x,a,1x,i8)') myid, trim(hostnam), idum
#endif

  close(unit=lunout)

#ifdef __MPI
  call MPI_FINALIZE(ierr);  call ckerr('MPI_FINALIZE')
#endif

  deallocate(iwrite)

CONTAINS
  subroutine ckerr(callnam)

    character*(*) :: callnam

    if(ierr.ne.0) then
       write(0,*) ' error ierr=',ierr,' in call "',trim(callnam), &
            '" on proc #',myid
       call bad_Exit
    endif

  end subroutine ckerr

end program proclist
