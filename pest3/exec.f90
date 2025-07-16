subroutine pstExec

  ! Execute method 

  USE pstcom
  IMPLICIT NONE
  integer,parameter :: r8 = SELECTED_REAL_KIND(12,100)

  CALL pstfunint

  ! find all rational surfaces
  CALL pstratsu2

  imode = 1
  do while(mm(imode) > 0 .and. imode <=nrun)

     if (imode>1) isolver=0

     ! reset small/large solutions & energies
     ! read profiles
     CALL pstreset

     ! radial mesh generation
     CALL pstnodes

     ! indexing setup
     CALL pstindint

     ! remap into new mesh
     CALL pstremap2

     ! finite element matrix assembly
     CALL pstmatelt
     
     ! add vacuum contribution
     CALL pstent34

     ! solve linear system
     CALL pstent23

     ! repeat any calculations ?
     CALL pstrepet2

     imode = imode + 1
  enddo

end subroutine pstExec
