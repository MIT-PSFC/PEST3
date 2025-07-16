#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function idescr_short(i)
      implicit none
      integer*2 i
      idescr_short = 0
      return
      end
 
      integer function idescr_long(i)
      implicit none
      integer i
      idescr_long = 0
      return
      end
 
      integer function idescr_float(f)
      implicit none
      real f
      idescr_float = 0
      return
      end
 
      integer function idescr_double(d)
      implicit none
      real*8 d
      idescr_double = 0
      return
      end
 
      integer function idescr_cstring(c)
      implicit none
      character c*(*)
      idescr_cstring = 0
      return
      end
 
      integer function idescr_longarr(i, idims, n)
      implicit none
      integer i
      integer idims(*), n
      idescr_longarr = 0
      return
      end
 
      integer function idescr_shortarr(i, idims, n)
      implicit none
      integer*2 i
      integer idims(*), n
      idescr_shortarr = 0
      return
      end
 
      integer function idescr_floatarr(f, idims, n)
      implicit none
      real f
      integer idims(*), n
      idescr_floatarr = 0
      return
      end
 
      integer function idescr_doublearr(d, idims, n)
      implicit none
      real*8 d
      integer idims(*), n
      idescr_doublearr = 0
      return
      end
 
      integer function idescr_cstringarr(c, idims, n)
      implicit none
      character *(*) c
      integer idims(*), n
      idescr_cstringarr = 0
      return
      end
#else
 
      integer function idescr_short(i)
      implicit none
      integer*2 i
      integer descr, mds__short
      idescr_short = descr(mds__short(), i, 0)
      return
      end
 
      integer function idescr_long(i)
      implicit none
      integer i
      integer descr, mds__long
      idescr_long = descr(mds__long(), i, 0)
      return
      end
 
      integer function idescr_float(f)
      implicit none
      real f
      integer descr, mds__float
      idescr_float = descr(mds__float(), f, 0)
      return
      end
 
      integer function idescr_double(d)
      implicit none
      real*8 d
      integer descr, mds__double
      idescr_double = descr(mds__double(), d, 0)
      return
      end
 
      integer function idescr_cstring(c)
      implicit none
      character c*(*)
      integer descr, mds__cstring, l
      idescr_cstring = descr(mds__cstring(), c, 0, len(c))
      return
      end
 
      integer function idescr_shortarr(i, idims, n)
      implicit none
      integer*2 i
      integer idims(*), n
      integer descr, mds__short
      integer j, jdims(8), num
      num = min(n,8)
      do j=1,num
         jdims(j) = idims(j)
      end do
      do j=num+1,8
         jdims(j) = 0
      end do
      idescr_shortarr = descr(mds__short(), i,
     $     jdims(1), jdims(2), jdims(3), jdims(4),
     $     jdims(5), jdims(6), jdims(7), jdims(8), 0)
      return
      end
 
      integer function idescr_longarr(i, idims, n)
      implicit none
      integer i
      integer idims(*), n
      integer descr, mds__long
      integer j, jdims(8), num
      num = min(n,8)
      do j=1,num
         jdims(j) = idims(j)
      end do
      do j=num+1,8
         jdims(j) = 0
      end do
      idescr_longarr = descr(mds__long(), i,
     $     jdims(1), jdims(2), jdims(3), jdims(4),
     $     jdims(5), jdims(6), jdims(7), jdims(8), 0)
      return
      end
 
      integer function idescr_floatarr(f, idims, n)
      implicit none
      real f
      integer descr, mds__float
      integer idims(*), n
      integer j, jdims(8), num
      num = min(n,8)
      do j=1,num
         jdims(j) = idims(j)
      end do
      do j=num+1,8
         jdims(j) = 0
      end do
      idescr_floatarr = descr(mds__float(), f,
     $     jdims(1), jdims(2), jdims(3), jdims(4),
     $     jdims(5), jdims(6), jdims(7), jdims(8), 0)
      return
      end
 
      integer function idescr_doublearr(d, idims, n)
      implicit none
      real*8 d
      integer descr, mds__double
      integer idims(*), n
      integer j, jdims(8), num
      num = min(n,8)
      do j=1,num
         jdims(j) = idims(j)
      end do
      do j=num+1,8
         jdims(j) = 0
      end do
      idescr_doublearr = descr(mds__double(), d,
     $     jdims(1), jdims(2), jdims(3), jdims(4),
     $     jdims(5), jdims(6), jdims(7), jdims(8), 0)
      return
      end
 
      integer function idescr_cstringarr(c, idims, n)
      implicit none
      character*(*) c
      integer descr, mds__cstring
      integer idims(*), n
      integer j, jdims(10), num
      num = min(n,8)
      do j=1,num
         jdims(j) = idims(j)
      end do
      jdims(num+1) = 0
      jdims(num+2) = len(c)
      idescr_cstringarr = descr(mds__cstring(), c,
     $     jdims(1), jdims(2), jdims(3), jdims(4),
     $     jdims(5), jdims(6), jdims(7), jdims(8),
     $     jdims(9), jdims(10))
      return
      end
#endif
 
 
 
 
