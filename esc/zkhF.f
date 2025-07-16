c-----------------------------------------------------      
c Fortran openning file
c-----------------------------------------------------      
      subroutine openftn(fname1)
      character fname1*1, fname*80
      integer i
      i=1
      do while(ichar(fname1(i:i)).ne.0 .and. i .lt. 14)
         i=i+1
      enddo
      fname=fname1(1:i)
      open(10,file=fname,form='unformatted',status='unknown')
      write(*,'(a,a)')'Output file=',fname
      end

c-----------------------------------------------------      
c Fortran closing file
c-----------------------------------------------------      
      subroutine closeftn
      close(10)
      end

      integer function ioftns(iU,fname,ch)
      character fname*80,ch
cc      if(ch.EQ.'r')then
cc      open(iU,file=fname,form='unformatted',recordtype='stream',iostat=i
cc     &o,status='old')
cc      else
cc      open(iU,file=fname,form='unformatted',recordtype='stream',iostat=i
cc     &o,status='unknown')
cc      endif
      if(ch.EQ.'r')then
         open(iU,file=fname,form='unformatted'
     &        ,iostat=io,status='old')
      else
         open(iU,file=fname,form='unformatted'
     &        ,iostat=io,status='unknown')
      endif
      if(io.NE.0)then
      ioftns=0
      else
      ioftns=1
      endif
      end
      integer function ioftnb(iU,fname,ch)
      character fname*80,ch
      if(ch.EQ.'r')then
      open(iU,file=fname,form='unformatted',iostat=io,status='old')
      else
      open(iU,file=fname,form='unformatted',iostat=io,status='unknown')
      endif
      if(io.NE.0)then
      ioftnb=0
      else
      ioftnb=1
      endif
      end
      integer function ioftna(iU,fname,ch)
      character fname*80,ch
      if(ch.EQ.'r')then
      open(iU,file=fname,iostat=io,status='old')
      else
      open(iU,file=fname,iostat=io,status='unknown')
      endif
      if(io.NE.0)then
      ioftna=0
      else
      ioftna=1
      endif
      end
      
      subroutine clftnb
      close(11)
      end
      subroutine clftna
      close(10)
      end
      
      subroutine fwri(i,n)
      integer i(n),n,k
      write(10)(i(k),k=1,n)
      end
      
      subroutine fwrd(x,n)
      double precision x(n)
      integer n,k
      write(10)(x(k),k=1,n)
      end
      
      subroutine fwrf(x,n)
      real*4 x(n)
      integer n,k
      write(10)(x(k),k=1,n)
      end
      
      subroutine frdi(i,n)
      integer i(n),n,k
      read(10)(i(k),k=1,n)
      end
      
      subroutine frdd(x,n)
      double precision x(n)
      integer n,k
      read(10)(x(k),k=1,n)
      end
      
      subroutine frdf(x,n)
      real*4 x(n)
      integer n,k
      read(10)(x(k),k=1,n)
      end

      subroutine openpp(i,q0)
      double precision q0
      real s
      character file1*23
      file1='/tmp/zakh2ahg/pp_01_1.0'
      write(file1(18:19),'(i2.2)')i
      s=q0
      write(file1(21:23),'(f3.1)')s
      open(10,file=file1,form='unformatted',status='unknown')
      write(*,'(a,a)')'open ',file1
      end
      subroutine esc_check
      character file1*23
      real*8 gy_i(25),q(25),F(25),Fs(25),Ps(25)
      file1='/tmp/zakh2ahg/pp_01_1.0'
1     write(*,'(a)',ADVANCE='NO')'Alan file N=>'
      read(*,*)i
      if(i.LT.0)stop
      write(file1(18:19),'(i2.2)')i
2     write(*,'(a)',ADVANCE='NO')'q_0=>'
      read(*,*)q0
      if(q0.LT.0)goto 1
      write(file1(21:23),'(f3.1)')q0
      open(10,file=file1,form='unformatted',status='unknown')
      write(*,'(a,a)')'open ',file1
      na1=25
      read(10)(gy_i(i),i=1,na1)
      read(10)(q(i),i=1,na1)
      read(10)(F(i),i=1,na1)
      read(10)(Fs(i),i=1,na1)
      read(10)(Ps(i),i=1,na1)
      close(10)
      open(2,file='/tmp/p1',status='unknown')
      da=1./(na1-1)
      do i=1,na1
      write(2,'(i4,1p6e12.4)')i,da*(i-1),gy_i(i),q(i),F(i),Fs(i),Ps(i)
      enddo
      close(2)
      goto 2
      end
      
      subroutine checkahg(r,z,psi,q,btxr,dbtxr,dp)
      character file1*79
      real*8 r(*),z(*),psi(*),q(*),btxr(*),dbtxr(*),dp(*)
      write(*,*)'file name'
      read(*,'(a)')file1
      open(10,file=file1,form='unformatted',status='unknown')
      write(*,'(a,a)')'open ',file1
      read(10)nr,ntau
      do ir=1,nr
      i=ntau*(ir-1)
      read(10)(r(itau+i),itau=1,ntau)
      read(10)(z(itau+i),itau=1,ntau)
      enddo
      read(10)(psi(ir),ir=1,nr)
      read(10)(q(ir),ir=1,nr)
      read(10)(btxr(ir),ir=1,nr)
      read(10)(dbtxr(ir),ir=1,nr)
      read(10)(dp(ir),ir=1,nr)
      close(10)
      end
C* :1 * 
      
