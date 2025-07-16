#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_value(expr, answer, retlen)
      implicit none
      character*(*) expr
      integer       answer, retlen
      mds_value = 0
      end
#else
c-----------------------------------------------------------------------
c mds_value
c get data
c
c Input
c   expr  -- string containing mds expression
c Output
c   answer -- answer descriptor
c   retlen -- return length
c Returns
c   1            -- success
c   even number  -- vms-like error
c
c 03/12/1999 ler Converting to unix usable
c 03/23/1999 ler Get return value by using size
c 04/12/1999 ler mds_value now handles return length
c-----------------------------------------------------------------------
      integer function mds_value(expr, answer, retlen)
      implicit none
      include 'MDSLIB'
 
      character*(*) expr
      integer       answer, retlen
 
      character*4097 str_expr
      integer status
 
      call mds__copywnull(str_expr, expr)
      status = MdsValue(str_expr, answer, 0, retlen)
 
      mds_value = status
 
      return
      end
#endif
 
 
 
 
 
 
