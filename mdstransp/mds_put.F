#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function mds_put(node, expr, dsc_data, units)
      implicit none
      character*(*) node
      character*(*) expr
      integer       dsc_data
      character*(*) units
      mds_put = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_put
c
c Input
c   node      -- character  i.e.
c   expr      -- character
c   dsc_data  -- descriptor pointing to array
c   units     -- character string
c
c Returns
c   vms style error code - odd numbers good
c                          even numbers bad
c   (odd numbers mapped to 1)
c
c **** WARNING *******
c   Is floating point convert stuff needed here?
c      For Digital Alpha? Digital VMS?, Solaris?
c ********************
c
c 03/17/1999 ler Converting to unix usable
c-----------------------------------------------------------------------
      integer function mds_put(node, expr, dsc_data, units)
      implicit none
 
      character*(*) node
      character*(*) expr
      integer       dsc_data
      character*(*) units
 
      character*4097 str_node, str_expr
      integer ln, le, lu
      integer dsc_units
      integer status
 
      integer descr, MdsPut
      integer mds__leniftrim, mds__cstring
 
      lu  = max(1, mds__leniftrim(units))
 
      call mds__copywnull(str_node, node)
      call mds__copywnull(str_expr, expr)
 
      dsc_units = descr(mds__cstring(), units, 0, lu)
      status = MdsPut(str_node, str_expr, dsc_data, dsc_units,0)
 
      mds_put = status
 
      return
      end
#endif
 
 
 
 
 
 
 
