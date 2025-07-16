C-------------------------------------------------------------------
      subroutine vprotec(tokstr,istat)
C
C  vprotec -- verify access permission
C
C    dmc Mar 1995 -- except at PPPL this routine should not do anything!
C      here at PPPL it is requested to use this to prevent general
C      TRANSP user access to JT60 data.
C
C    rewritten DMC Nov 1996 -- access list now in a file.
C      assume wide open if cannot find file.
C
      implicit NONE
C
C-------------------------------------
C
      character*(*) tokstr  ! tok or tok.yy string, input
      integer istat         ! status code, output. istat=0:  OK
C
C-------------------------------------
C
      integer numacl,icln,il,itok,ilnb
C
      character*80 zbuff,zline
C
      character*4 tokacl(100)
      logical     grantd(100)
C
      integer str_length,ilen
C
      data numacl/-1/
C
      save tokacl
      save grantd
C
C-------------------------------------
C
      istat=0
C
C  check for existence of ACL file; collect list of ACL protected
C  tokamak names
#ifdef __VMS
      if(numacl.eq.0) go to 9999
      if(numacl.eq.-1) then
         numacl=0
         open(unit=99,file='TRANSP$:[USER]RPLOT.ACL',readonly,
     >        status='old',err=9999)
 10      continue
         read(99,'(A)',end=19) zline
         if(zline(1:1).eq.'*') then
            icln=index(zline,':')
            numacl=numacl+1
            tokacl(numacl)=' '
            tokacl(numacl)=zline(2:icln-1)
            grantd(numacl)=.false.
         endif
         go to 10
 19      close(unit=99)
      endif
C
C  extract passed tokamak name
C
      il=index(tokstr,'.')-1
      if(il.lt.0) il=index(tokstr,']')-1
      if(il.lt.0) il=index(tokstr,'>')-1
      if(il.lt.0) il=index(tokstr,' ')-1
      if(il.lt.0) il=len(tokstr)
C
C  is it on the list of tokamaks?
C
      do itok=1,numacl
         if(tokstr(1:il).eq.tokacl(itok)) go to 20
      enddo
      go to 9999
C
 20   continue
      if(grantd(itok)) go to 9999
      write(6,
     >  '('' %vprotec:  checking access list for '',a)') tokacl(itok)
C
C  get username
C
      ilen=0
      zbuff=' '
      call sget_user(zbuff)
      ilen=str_length(zbuff)
C
C  restrict access to users named here
C
      open(unit=99,file='TRANSP$:[USER]RPLOT.ACL',readonly,
     >     status='old',err=9999)
C
 100  continue
      read(99,'(A)',end=119) zline
      if(zline(1:1).ne.'*') go to 100
      icln=index(zline,':')
      if(zline(2:icln-1).ne.tokacl(itok)) go to 100
C
 110  continue
      read(99,'(A)',end=119) zline
      if(zline(1:1).eq.'*') go to 119
      ilnb=max(1,(index(zline,' ')-1))
      if(zline(1:ilnb).eq.zbuff(1:ilen)) go to 120
      go to 110
C
 119  continue
      close(unit=99)
      istat=1
      write(6,9901) tokacl(itok),zbuff(1:ilen)
      go to 9999
C
 120  continue
      grantd(itok)=.TRUE.
      close(unit=99)
C
C------------------------
 9901          format(/' ?access to ',a,' data restricted'/
     >                 '  user "',a,'" not authorized.'/
     >                 '  Email transp@pppl.gov for help')
C
 9999 continue
#endif
      return
      end
