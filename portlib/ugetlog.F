      subroutine ugetlog(zbuff)
      character*(*) zbuff               ! for username string, output
c
c  return the current username
c
#ifdef __UNIX
c
      character*80 homedir
      character*150 ztmpfile
      integer ilz,iln,str_length,ier
c
      call sget_env('LOGNAME',zbuff)
      if(zbuff.eq.' ') then
         call sget_env('USER',zbuff)
         if(zbuff.eq.' ') then
c
c  (perhaps the C routine getlogin_r() could be used; a portlib C wrapper
c  would have to be added; BYTE_DECLARE used, etc.)
c
            write(0,*) ' ?ugetlog,sget_user: could not get username.'
            zbuff = 'unknown'

         endif
      endif
#endif
      return
      end
