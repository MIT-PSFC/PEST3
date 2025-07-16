      subroutine get_arg_count(numargs)
C
#if __NAGWARE
      use f90_unix_env
#endif
C
      integer numargs
C
C  return the number of arguments on the command line
C
#ifdef __VMS
      numargs=iargc_vms(numargs)
#endif
 
#ifdef __UNIX
 
#if !__CRAY && !__NAGWARE
      integer iargc
#elif __CRAY || __X1
      external ipxfargc
      integer ipxfargc
#endif
 
#if __CRAY || __X1
      numargs = ipxfargc()
#elif __ABS
      numargs = iargc_()
#else
      numargs = iargc()
#endif
 
#endif /* UNIX */
 
      return
      end
 
#ifdef __LF95

c
c  this compiler lacks FY-2003 routine -- attempt a patch...
c

      integer function command_argument_count()

      integer :: iargs

      call get_arg_count(iargs)
      command_argument_count = iargs

      return
      end

#endif
