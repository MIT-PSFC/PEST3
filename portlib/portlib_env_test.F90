program portlib_env_test

  ! test portlib modules to support MPI distribution of working directory
  ! and command line arguments and selected environment variables

  ! May 19, 2011: added testing of MPI proc/node information (mpi_proc_data)
  ! which for MPI_COMM_WORLD is set up on the first mpi_set_env or mpi_get_env
  ! call or by mpi_share_env

  use logmod
  use mpi_proc_data
  implicit none

#ifdef __MPI
  include 'mpif.h'
#endif

  ! local communicator
#ifdef __MPI
  integer :: nbi_comm = MPI_COMM_WORLD
#else
  integer :: nbi_comm = -1
#endif

  integer :: myid,ierr,iertmp,istat,jstat,iolun,idot,iost,ierloc,io
  integer :: num_nodes,maxlist,mylist_size,itest,mynode
  integer, dimension(:), allocatable :: myprocs,myproc2

  character*200 :: env_val,log_file,jobid
  
  character*32 runid
  character*200 :: log_file_prefix,zcmd,path_transl

  character*40 filexxx,filennn,zpbsid

  logical :: pflag,pexist
  integer :: jsystem,izp,izleft,iotest

#ifdef __MPI
  character*(MPI_MAX_PROCESSOR_NAME) proc_name
  integer iresult
#endif
  integer :: numprocs = 1

  !-----------------------------------------
  ! executable code

#ifdef __MPI

  call MPI_INIT(ierr)
  if(ierr.ne.0) call errmsg_exit(' MPI_INIT failure! ')

  call MPI_COMM_RANK(nbi_comm,myid,ierr);   call ckerr('MPI_COMM_RANK')
  call MPI_COMM_SIZE(nbi_comm,numprocs,ierr);   call ckerr('MPI_COMM_SIZE')
  if(myid.eq.0) then
     write(0,*) ' mpi portlib_env_test:  numprocs = ',numprocs
     write(0,*) ' communicator ID:       nbi_comm = ',nbi_comm
  endif

  call MPI_GET_PROCESSOR_NAME(proc_name,iresult,ierr)
  call ckerr('MPI_GET_PROCESSOR_NAME')

  write(0,*) ' portlib_env_test: MPI, myid=',myid,' proc_name = "', &
       trim(proc_name),'"'
#else
  myid = 0
  write(0,*) ' portlib_env_test: serial operation (myid=0).'
#endif

  if(myid.eq.0) then
     call mpi_sget_env('PBS_JOBID',env_val,ierr); call ckerr('mpi_sget_env(1)')
     idot = index(env_val,".") - 1
     if(idot.gt.1) then
        call mpi_sset_env('PBS_JOBID',env_val(1:idot),ierr)
        call ckerr('mpi_sset_env(1)')
     endif
  endif
  call mpi_share_env(myid,ierr); call ckerr('mpi_share_env')
  if(myid.eq.0) then
     call mpi_printenv(myid,6)
  else
     
     call mpi_sget_env('PBS_JOBID',jobid,ierr); call ckerr('mpi_sget_env')
     
  endif
  log_file = "env_test.log"
  call openLog(file=trim(log_file),status='ONCE',iostat=iost, &
       comm_world=nbi_comm) 

  call setLogLevel('info') ! log file for info - error messages

  call infoLog('portlib_env_test log file opening')
  call infoLog('portlib_env_test log file opening')
  call infoLog('portlib_env_test log file opening')

  if(myid.eq.0) then

     ! uncomment this for extra debugging in mpi_share_env(...):
     !!! call mpi_sset_env('MPI_SHARE_ENV_DEBUG','TRUE',ierr)
     call mpi_sget_env('PBS_JOBID',env_val,ierr); call ckerr('mpi_sget_env(1)')
     idot = index(env_val,".") - 1
     if(idot.gt.1) then
        call mpi_sset_env('PBS_JOBID',env_val(1:idot),ierr)
        call ckerr('mpi_sset_env(1)')
     endif
     call mpi_sget_env('PREACTDIR',env_val,ierr); call ckerr('mpi_sget_env(2)')
     call mpi_sset_env('NUBEAM_ACTION','test',ierr); call ckerr('mpi_sset_env(2)')
     call mpi_sset_env('NUBEAM_WORKPATH','foo',ierr); call ckerr('mpi_sset_env(3)')

  endif

  call mpi_share_env(myid,ierr); call ckerr('mpi_share_env (1)')

  if(myid.eq.0) then
     call mpi_printenv(myid,6)   ! output the environment 
     call mpi_sset_env('NUBEAM_COMPONENT','TRUE',ierr); call ckerr('mpi_sset_env(4)')
     call mpi_sset_env('TORIC_COMPONENT','FALSE',ierr); call ckerr('mpi_sset_env(4)')
     call mpi_sset_env('PT_SOLVER_COMPONENT','FALSE',ierr); call ckerr('mpi_sset_env(4)')
  endif

  call mpi_share_env(myid,ierr); call ckerr('mpi_share_env (2)')

  if(myid.gt.0) then
     iolun=6  ! 6 is OK -- was tested
  else
     iolun=6
  endif

#ifdef __MPI
  if(myid.le.0) then
     write(iolun,*) &
          ' child logs written to: env_test_<jobid>.log_1, env_test_<jobid>.log_2, ... '
     write(iolun,*) '   in directory matching cwd of parent.'
  endif
#endif

  if(myid.gt.0) then

     call mpi_sget_env('PBS_JOBID',jobid,ierr); call ckerr('mpi_sget_env')
     log_file = ' '
     if(myid.le.9) then
        write(log_file,'(a,i1)') "env_test_"//trim(jobid)//".log_",myid
     else if(myid.le.99) then
        write(log_file,'(a,i2)') "env_test_"//trim(jobid)//".log_",myid
     else if(myid.le.999) then
        write(log_file,'(a,i3)') "env_test_"//trim(jobid)//".log_",myid
     else if(myid.le.9999) then
        write(log_file,'(a,i4)') "env_test_"//trim(jobid)//".log_",myid
     endif

     open(unit=iolun,file=trim(log_file),status='unknown',iostat=istat)
     if(istat.eq.0) then
        write(0,*) ' proc #',myid,' opened file OK: ',trim(log_file)
     else
        write(0,*) ' ?? proc #',myid,' file open failure: ',trim(log_file)
     endif
  endif

  call mpi_printenv(myid,iolun)   ! output the environment 

  if(myid.gt.0) then
     close(unit=iolun,iostat=istat)
     if(istat.eq.0) then
        write(0,*) ' proc #',myid,' closed file OK: ',trim(log_file)
     else
        write(0,*) ' ?? proc #',myid,' file close failure: ',trim(log_file)
     endif
  endif

#ifdef __MPI
  call MPI_BARRIER(nbi_comm,ierloc)
#endif

  call mpi_print_info(nbi_comm,0,ierloc)
  if(numprocs.gt.1) call mpi_print_info(nbi_comm,1,ierloc)
  if(numprocs.gt.2) call mpi_print_info(nbi_comm,numprocs-1,ierloc)

  itest = max(0,numprocs-1)

  if(myid.eq.itest) then
     call mpi_proc_info(nbi_comm,ierloc, num_nodes=num_nodes, &
          maxnode_numprocs=maxlist,mynode=mynode)

     allocate(myprocs(maxlist+1)) ! deliberately padded
     allocate(myproc2(maxlist+1)) ! deliberately padded

     call mpi_proc_info(nbi_comm,ierloc, node_proclist=myprocs, &
          node_numprocs = mylist_size)

     call mpi_proc_info(nbi_comm,ierloc, node_number=mynode, &
          node_proclist=myproc2, &
          parallel_file_system = pflag)

     write(0,*) ' '
     write(0,*) ' these should match: '
     write(0,*) '    myprocs: ',myprocs(1:mylist_size+1)
     write(0,*) '    myproc2: ',myproc2(1:mylist_size+1)
     write(0,*) ' '
     write(0,*) ' parallel_file_system: ',pflag
     write(0,*) ' '
  endif

  call closeLog

#ifdef __MPI
  call MPI_BARRIER(nbi_comm,ierr)
  if(ierr.ne.0) write(0,*) ' ?? proc#',myid,' MPI_BARRIER returned ierr=',ierr
#endif

  ! test file copying
  filexxx = 'testxxxxx.dat'
  call find_io_unit(io)
  if(myid.eq.0) then
     open(unit=io,file=filexxx,status='unknown')
     write(io,*) 'this is...'
     write(io,*) '  ...a test'
     close(unit=io)
  endif

  do io=0,numprocs-1
     filennn=' '
     write(filennn,'("test",i5.5,".dat")') io
     call mpi_file_bcast(nbi_comm,filexxx,ierr, newname=filennn)
     if(ierr.ne.0) write(0,*) ' ?? proc#',myid,' parcopy returned ierr=',ierr
  enddo

#ifdef __MPI
  !! write(0,*) ' enter final MPI_ALLREDUCE: myid = ',myid
  call MPI_ALLREDUCE(ierr,istat,1,MPI_INTEGER,MPI_MAX,nbi_comm,iertmp)
  !! write(0,*) ' after final MPI_ALLREDUCE: myid=',myid,' iertmp=',iertmp
  write(0,*) ' at MPI_FINALIZE: myid = ',myid,' ierr = ',ierr
  call MPI_FINALIZE(ierr);  call ckerr('MPI_FINALIZE')
#endif

  if(myid.eq.0) then

     istat = jsystem('mkdir -p foonzork')
     write(0,*) 'jsystem("mkdir -p foonzork") status: ',istat
     istat = jsystem('rmdir foonzork')
     write(0,*) 'jsystem("rmdir foonzork") status: ',istat

     inquire(file='/usr/pppl/bin/pbs_timeleft',exist=pexist,iostat=istat)
     pexist = pexist .AND. (istat.eq.0)
     if(.NOT.pexist) then
        write(0,*) ' pbs_timeleft command script not found: pexist=',pexist, &
             ' istat= ',istat
     else
        call sget_env('PBS_JOBID',zpbsid)
        if(zpbsid.eq.' ') then
           write(0,*) ' pbs_timeleft not attempted: PBS_JOBID not found.'
        else
           write(0,*) ' PBS_JOBID: '//trim(zpbsid)
           izp=index(zpbsid,'.')-1
           if(izp.le.0) izp=len_trim(zpbsid)

           zcmd = '/usr/pppl/bin/pbs_timeleft '//zpbsid(1:izp) &
                //' > loc_'//zpbsid(1:izp)//'.tmp'

           write(0,*) ' attempting: jsystem("'//trim(zcmd)//'")'

           istat = jsystem(zcmd)

           if(istat.ne.0) then
              write(0,*) ' jsystem failed, error status: ',istat
           else
              !  read pbs_timeleft output
              call find_io_unit(iotest)
              open(unit=iotest,file='loc_'//zpbsid(1:izp)//'.tmp', &
                   status='old',iostat=jstat)
              if(jstat.ne.0) then
                 write(0,*) ' jsystem OK but no file: loc_'//zpbsid(1:izp)//'.tmp'
              else
                 read(iotest,*,iostat=jstat) izleft
                 if(jstat.ne.0) then
                    write(0,*) ' error reading:  loc_'//zpbsid(1:izp)//'.tmp'
                    close(iotest)
                 else
                    close(iotest,status='delete')
                    write(0,*) ' pbs_timeleft OK, seconds=',izleft, &
                         ' hours=',izleft/3600
                 endif  ! pbs_timeleft output file read OK
              endif  ! pbs_timeleft output file open OK
           endif  ! jsystem OK
        endif  ! PBS_JOBID defined
     endif  ! PPPL pbs_timeleft procedure found

     write(0,*) ' '
     call nodepath_translator('/u/dmccune',path_transl,nbi_comm,0)
     write(0,*) ' /u/dmccune path translation: '//trim(path_transl)
     call nodepath_translator('/local/dmccune',path_transl,nbi_comm,0)
     write(0,*) ' /local/dmccune path translation: '//trim(path_transl)

  endif  ! myid.eq.0

CONTAINS
  subroutine ckerr(callnam)

    character*(*) :: callnam

    if(ierr.ne.0) then
       write(0,*) ' error ierr=',ierr,' in call "',trim(callnam), &
            '" on proc #',myid
       call bad_Exit
    endif

  end subroutine ckerr

end program portlib_env_test
