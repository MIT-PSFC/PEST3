#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_errstr(ierr, str)
      implicit none
      integer ierr
      character*(*) str
      str = ' '
      mds_errstr = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_errstr
c get error string
c
c 03/16/1999 ler Converting to unix usable
c                Get rid of null end of string before returning
c-----------------------------------------------------------------------
        integer function mds_errstr(ierr, str)
        implicit none
 
        integer ierr
        character*(*) str
	
        integer dsc, size, status
        character*19 cmd
 
        integer descr, MdsValue, mds__cstring
 
        str=' '
        write(cmd, '(''GETMSG('',I11,'')'')') ierr
        dsc = descr(mds__cstring(), str, 0, len(str))
        status = MdsValue(cmd //char(0), dsc, 0, size)
        call mds_addblanks(str)
 
        mds_errstr = status
 
        return
        end
#endif
 
 
