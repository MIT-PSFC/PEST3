      subroutine eqm_premsg(zmsg)
C
      use eq_module
C
C  set the pre-message for xplasma error handling routines
C  if this is set, each error message is preceded with this message text.
C
      IMPLICIT NONE
C
      character*(*) zmsg
C
c$cface_input:  zmsg
C------------------------------------
C
      eq_premsg_text=zmsg
C
      return
      end
