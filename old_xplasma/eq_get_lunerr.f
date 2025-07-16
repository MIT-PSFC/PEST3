      subroutine eq_get_lunerr(ilunerr)
!
!  return the lun used for eq/xplasma error messages
!
      use eq_module
!
      ilunerr = lunerr
!
      return
      end
