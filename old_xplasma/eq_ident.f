      subroutine eq_ident(eqmodel_o,axisymm_o,scrapeoff_o)
C
      use xplasma_obj_instance
      use eq_module
C
C  return identification of equilibrium data model as currently stored
C
      IMPLICIT NONE
C
      character*(*) eqmodel_o           ! name of model
      integer axisymm_o                 ! axisymmetry flag, =1 for axisymmetry
      integer scrapeoff_o               ! scrapeoff flag, =1: scrapeoff layer
C
C-----------------------------------
      integer :: ier
      logical :: flag
C-----------------------------------
C
      call xplasma_global_info(s,ier, initLabel=eqmodel_o)
      if(ier.ne.0) then
         write(lunerr,*) ' %eq_ident error (xplasma) ier=',ier
         call xplasma_error(s,ier,lunerr)
         axisymm_o=1
         scrapeoff_o=0
      else
         call xplasma_global_info(s,ier,axisymm=flag)
         if(flag) then
            axisymm_o=1
         else
            axisymm_o=0
         endif
         call xplasma_global_info(s,ier,scrapeoff=flag)
         if(flag) then
            scrapeoff_o=1
         else
            scrapeoff_o=0
         endif
      endif
C
      return
      end
