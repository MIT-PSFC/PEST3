      subroutine vmsieee(lun_warn,zvals,ilim,ufieee,ieee)
c
c  VMS:  convert single precision floating point data to IEEE format
c  UNIX:  just copy the numbers
c
c  overflow/underflow check provided.
c
      integer lun_warn                  ! LUN for warning messages
      integer ilim
      real zvals(ilim)                  ! data to convert or recover
      real ufieee(ilim)                 ! converted data (ieee format)
      integer ieee                      ! =1:  convert to ieee
c                                       ! =0:  convert from ieee back
#ifdef __VMS
      INTEGER CVT$K_VAX_F,CVT$K_IEEE_S
      PARAMETER (CVT$K_VAX_F=0,CVT$K_IEEE_S=4)
#endif
c
      iwarnuf=0
      iwarnof=0
c
      do i=1,ilim
         if(ieee.eq.1) then
C
C  convert to ieee
C  check for over/underflow
C
            if(zvals(i).ne.0.0) then
               if(abs(zvals(i)).lt.1.1754945e-38) then
                  zvals(i)=0.0
                  if((iwarnuf.eq.0).and.(lun_warn.gt.0)) then
                     write(lun_warn,
     >                  '('' %ubwfcmp:  underflow warning.'')')
                     iwarnuf=1
                  endif
               else if(zvals(i).gt.8.4070587e+37) then
                  zvals(i)=8.4070587e+37
                  if((iwarnof.eq.0).and.(lun_warn.gt.0)) then
                     write(lun_warn,
     >                  '('' %ubwfcmp: +overflow warning.'')')
                     iwarnof=1
                  endif
               else if(zvals(i).lt.-8.4070587e+37) then
                  zvals(i)=-8.4070587e+37
                  if((iwarnof.eq.0).and.(lun_warn.gt.0)) then
                     write(lun_warn,
     >                  '('' %ubwfcmp: -overflow warning.'')')
                     iwarnof=1
                  endif
               endif
            endif
#ifdef __UNIX
            ufieee(i)=zvals(i)
#endif
#ifdef __VMS
            call cvt$convert_float(zvals(i),%VAL(CVT$K_VAX_F),
     >         ufieee(i),%VAL(CVT$K_IEEE_S),%VAL(0))
#endif
         else
c
c  convert back *from* ieee; no error check
c
#ifdef __UNIX
            zvals(i)=ufieee(i)
#endif
#ifdef __VMS
            call cvt$convert_float(ufieee(i),%VAL(CVT$K_IEEE_S),
     >                 zvals(i),%VAL(CVT$K_VAX_F),%VAL(0))
#endif
c
         endif
      enddo
c
      return
      end
 
