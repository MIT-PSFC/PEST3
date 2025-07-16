subroutine execsystem_test
#ifdef __UNIX
  use execsystem

  integer      :: ier
  character*64 :: arg               ! single argument
  character*32 :: env0(0), env3(5)  ! new environments
  character*64 :: marg(6)           ! multiple arguments
  character*16 :: exec              ! executable

  ! -- initialize --
  EXEC_VERBOSE = -1   ! verbose messages to stderr

  env3(1) = "  DOG=RUN_AWAY"
  env3(2) = " NO_VALUE=  "
  env3(3) = "   CAT=mouse/meow"
  env3(4) = "=NO_KEY "
  env3(5) = "MOUSE=/time/to/roar"
     
  ! -- test --
  EXEC_STATIC  = 1
  arg = "  echo 'hello world with static arguments'"
  ier = jsystem(arg)
  print '(a,i5)', '%ptest: hello result = ',ier
  if (ier/=0) call bad_exit
   
  EXEC_STATIC  = 0
  arg = "  echo 'hello world with heap arguments'"
  ier = jsystem(arg)
  print '(a,i5)', '%ptest: hello result = ',ier
  if (ier/=0) call bad_exit
   
  exec    = "ls"
  marg(1) = "ls"
  marg(2) = "-rt"

  ier = jexec(exec, 2, marg)
  print '(a,i5)', '%ptest: ls result = ',ier
  if (ier/=0) call bad_exit

  exec    = "printenv"
  marg(1) = "printenv"

  ier = jexec(exec,1,marg,size(env3),env3)
  print '(a,i5)', '%ptest: exec printenv result = ',ier
  if (ier/=0) call bad_exit

  EXEC_VERBOSE = 0   ! turn off verbose messages
#endif
end subroutine execsystem_test
