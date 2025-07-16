#ifdef __NOMDSPLUS
#define STUB
#endif
#ifdef STUB
      integer function MDS_PUT3D(node, expr, dsc_data, units,
     $     dsc_dim0, udim0, dsc_dim1, udim1, dsc_dim2, udim2)
      implicit none
      character*(*) node
      character*(*) expr
      integer       dsc_data
      character*(*) units
      integer       dsc_dim0
      character*(*) udim0
      integer       dsc_dim1
      character*(*) udim1
      integer       dsc_dim2
      character*(*) udim2
      mds_put3d = 0
      return
      end
#else
c-----------------------------------------------------------------------
c mds_put3d
c
c Input
c   node      -- character  i.e.
c   expr      -- character
c   dsc_data  -- descriptor pointing to array
c   units     -- character string
c   dsc_dim0  -- descriptor pointing to array
c   udim0     -- character string
c   dsc_dim1  -- descriptor pointing to array
c   udim1     -- character string
c   dsc_dim2  -- descriptor pointing to array
c   udim2     -- character string
c
c Returns
c   vms style error code - odd numbers good
c                          even numbers bad
c   (odd numbers mapped to one)
c
c **** WARNING *******
c   Is floating point convert stuff needed here?
c      For Digital Alpha? Digital VMS?, Solaris?
c ********************
c
c 03/17/1999 ler Converting to unix usable
c-----------------------------------------------------------------------
      integer function mds_put3d(node, expr, dsc_data, units,
     $     dsc_dim0, udim0, dsc_dim1, udim1, dsc_dim2, udim2)
      implicit none
 
      character*(*) node
      character*(*) expr
      integer       dsc_data
      character*(*) units
      integer       dsc_dim0
      character*(*) udim0
      integer       dsc_dim1
      character*(*) udim1
      integer       dsc_dim2
      character*(*) udim2
 
      character*4097 str_node, str_expr, str_units
      character*4097 str_udim0, str_udim1, str_udim2
      integer ln, le, lu, lu0, lu1, lu2
      integer dsc_units, dsc_udim0, dsc_udim1, dsc_udim2
      integer status
 
      integer mds__leniftrim, mds__cstring
      integer descr, MdsPut
 
      ln  = max(1, mds__leniftrim(node))
      le  = max(1, mds__leniftrim(expr))
      lu  = max(1, mds__leniftrim(units))
      lu0 = max(1, mds__leniftrim(udim0))
      lu1 = max(1, mds__leniftrim(udim1))
      lu2 = max(1, mds__leniftrim(udim2))
 
      call mds__copywnull(node, str_node)
      call mds__copywnull(expr, str_expr)
      call mds__copywnull(units, str_units)
      call mds__copywnull(udim0, str_udim0)
      call mds__copywnull(udim1, str_udim1)
      call mds__copywnull(udim2, str_udim2)
 
      dsc_units = descr(mds__cstring(), str_units, 0, lu)
      dsc_udim0 = descr(mds__cstring(), str_udim0, 0, lu0)
      dsc_udim1 = descr(mds__cstring(), str_udim1, 0, lu1)
      dsc_udim2 = descr(mds__cstring(), str_udim2, 0, lu2)
      status = MdsPut(str_node, str_expr,
     $     dsc_data, dsc_units, dsc_dim0, dsc_udim0,
     $     dsc_dim1, dsc_udim1, dsc_dim2, dsc_udim2, 0)
 
      mds_put3d = status
 
      return
      end
#endif
