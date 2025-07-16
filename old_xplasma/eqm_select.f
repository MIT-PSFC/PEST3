      subroutine eqm_select(model,ilaxi)
C
      use xplasma_obj_instance
      use eq_module
C
      IMPLICIT NONE
C
C  select an equilibrium model type
C
      character*(*) model       ! input label string
      integer ilaxi             ! input =1 if axisymmetric, o/w FALSE
      integer :: ier
C
C----------------------------------
C
      if(.not.eq_module_init) then
         call xoi_init(ier)
         if(ier.ne.0) write(lunerr,*)
     >        ' ?eqm_select: xplasma_init returned ier = ',ier
         eq_module_init = .TRUE.
      endif
C
      call xplasma_init(s,model,ier)
      if(ier.ne.0) write(lunerr,*)
     >     ' ?eqm_select: xplasma_init returned ier = ',ier
C
      if(ilaxi.eq.1) then
         continue               ! axisymmetry is supported
      else
         write(lunerr,*) 
     >        ' ?eqm_select: non-axisymmetric XPLASMA not implemented!'
      endif
C
      return
      end
