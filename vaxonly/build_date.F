      subroutine build_date( cvsn_date, istat )
*.. build_date                          return transp common build date, 
*                                       from $TRANSP_BUILD_DATE, if set,
*                                       else from source/misc/build_date
*
*
*     Author                            jim.conboy@jet.uk
*     Version                           1.00,  26May2006
*     Modifications
*
*     1.00  26May2006                   jim.conboy@jet.uk
*                                       Extracted from mdsplot.for
*----_^---------------------------------=========================================|
      implicit  none
*
*
*..Arguments
*
      character(len=24),intent(out)  :: cvsn_date      ! transp common build date
      integer, intent(out)           :: istat          ! return status, 0=OK
*
*..local
      integer                        :: i,             ! loop index 
     +                                  iop,           ! file open opt
     +                                  luntrm=6       ! IO unit
*
      character(len=256)             :: zfile          ! file name buffer
*
      character(len=*),parameter     :: zcr='build_date'
*----_^------------------------------::================!.........................|
*
C
C  get TRANSP code vsn date string - may be provided in environment,
C                                    if we are lucky
C
      call sget_env( 'TRANSP_BUILD_DATE', cvsn_date )
      if(len_trim(cvsn_date) > 20 )                     then
         istat  =  0 ;                                 return
                                                       endif
*
*--    never mind, do it the hard way..
*
#ifdef __VMS
      i=1
      call ufopen(99,'source_misc:xshare_build.date',1,istat,0)
#else
      write(luntrm,*) zcr,': call sget_env'
      call sget_env('CODESYSDIR',zfile)
      i=index(zfile,' ')-1
      write(luntrm,*) zcr,': call ufopen xshare_build.dat'
      iop = 1
      call ufopen(99,zfile(1:i)//'/source/misc/xshare_build.date',
     >     iop, istat, 0)
      if(istat.ne.0) call ufopen(99,'xshare_build.date',iop,istat,0)
#endif
      if(istat.ne.0) then
         if (i .gt. 1) then
            write(luntrm,*)
     1           'open failure:  source_misc:xshare_build.date'
            call bad_exit
         else
            write(luntrm,*) zcr,':  build date unknown.'
            cvsn_date = 'unknown'
         endif
      else
         read(99,'(A)') cvsn_date
         close(unit=99)
      endif
*
      end subroutine build_date
