  subroutine logmod_test( cFile, cLvl, istat )
  !.. logmod_test                      test logmod
  !                                    codesys/source/porttest/porttest.logmod_test.
  !
  !
  !  Author                            jim.conboy@ccfe.ak.uk
  !  Version                           1.00,  25Jan2010
  !  Modifications
  !
  !  1.00   25Jan2010                  jim.conboy@ccfe.ak.uk
  !                                    From comments in logmod.f90
  !----_^---------------------------------=========================================|
  !
  use logmod
  !
  implicit  none
  !
  !..GLOBAL
  !variable                          :: name           ! description
  !
  !..Arguments
  !
  character(len=*),intent(in)        :: cFile          ! Name of File to open, else 'None' to skip
  character(len=*),intent(in)        :: cLvl           ! Info level, or 'df' for default
  integer, intent(out)               :: istat          ! return status, 0=OK
  !
  !..local
  integer                            :: i = 47
  real*8                             :: d = 3.14159e17
  character(len=32)                  :: cBuf           ! value of $LOG_LEVEL
  !
  character(len=*),parameter         :: cSr='logmod_test'
  !----------------------------------::================!.........................|
  !
  write(*,'(3a)') cSr,' -I- testing logmod, file ',trim(cFile),' level ',trim(cLvl)
  !
  if (isOpenLog()) then
     write(*,'(a)') '!logmod_test: The log file is OPEN'
  else
     write(*,'(a)') '!logmod_test: The log file is CLOSED'
  end if

  cBuf   = cFile(:4)
  call ulower(cBuf)
  if( cBuf(:4) /= 'none' )  &
     call openLog(cFile,status='REPLACE',icpu=7,numcpu=16)     ! if not opened, logging goes to portlib_<pid>.log
  !
  cBuf = cLvl
  call ulower(cBuf) 
  if( cBuf(:2) == 'df' )                          then
     call sget_env( 'LOG_LEVEL', cBuf    )
     write(*,'(3a)') cSr,' -I- Env. Var $LOG_LEVEL = ',cBuf 
  else
     call setLogLevel(cLvl)   ! default log level = WARN_LEVEL to show warnings and errors only
  endif
  !
  write(6,'(2a,i2)')  cSr,' -I- File Log level is ', loglevel()
  !
  call infoLog("just an info message")
  call infoLog("just an info message with i = ",i)
  call infoLog("just an info message with d = ",d)
  !
  call enterLog("logmod_test")    ! convenient place to put this
  !
  call infolog( cSr, 'Info with routine name ' )
  call infolog( cSr, 'Info with routine name, i = ',i )
  !
  call warnLog("just a warn message")
  call enterLog("logmod_test::short_section")   
  call gosleep(1)
  call warnLog("just a warn message with i = ",i)
  call warnLog("just a warn message with d = ",d)
  !
  call errorLog("just an error message")
  call errorLog("just an error message with i = ",i)
  call gosleep(1)
  call errorLog("just an error message with d = ",d)
  call exitLog("logmod_test::short_section")   
  !
  call exitLog("logmod_test")    ! convenient place to put this
  !   
  call warnLog("enable buffering and no time stamps")
  !
  if (logLevel()>=INFO_LEVEL) then
    call setLogBuffered(.true.)     ! do not flush after each log message
    call setLogTime(.false.)        ! do not time stamp log messages
  !
    call infoLog("some more detailed messages")
  !
    write(logline,'(a,i2,a,e12.4)') 'vec(',i,') = ',d
    call infoLog(logline)
  !
    call setLogTime(.true.)
    call setLogBuffered(.false.)     ! flush buffers
  end if
  !   
  call warnLog("time stamps enabled")
  !
  end subroutine logmod_test
