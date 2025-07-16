module xplasma_definitions

  ! public interface module for xplasma library
  ! kernel is in xplasma_obj; auxilliary routine interfaces defined here...

  use xplasma_fromgeq    ! xplasma_fromgeqdsk (separated for debugging)

  use xplasma_avg1dmod   ! 1d volume avg (separated for debugging)

  use xplasma_profs   ! routines that aid in creation of profile functions

  use xplasma_mcgrid  ! routines for MCgrid irregular grids & assoc. profiles.

  use xplasma_nclass  ! calculate Pfirsch-Schlutter moments for NCLASS

  ! also used by xplasma_profs:
  use xplasma_rzgeo   ! RZ flux surface geometry definition

  ! also used by xplasma_profs, xplasma_rzgeo, xplasma_mcgrid, xplasma_nclass:
  use xplasma_flxint  ! numerical integrator routines

  ! also used by xplasma_profs, xplasma_rzgeo:
  use xplasma_sol     ! scrape off layer data -- e.g. limiters

  ! also used by xplasma_profs, xplasma_sol, xplasma_rzgeo:
  use xplasma_ctran   ! coordinate transformations & extrema

  !-------------------------------------
  !  **kernel** used by all
  use xplasma_obj     ! xplasma data types & kernel routines

end module xplasma_definitions
