      SUBROUTINE scrfree
      USE mintrp
      USE scrunch_inc1
      USE scrunch_inc2
!...begin restore memory
       IF( ALLOCATED(rin) ) DEALLOCATE(rin)
       IF( ALLOCATED(zin) ) DEALLOCATE(zin)
       IF( ALLOCATED(xvec) ) DEALLOCATE(xvec)
       IF( ALLOCATED(gvec) ) DEALLOCATE(gvec)
       IF( ALLOCATED(xdot) ) DEALLOCATE(xdot)
       IF( ALLOCATED(xstore) ) DEALLOCATE(xstore)
       IF( ALLOCATED(m1) ) DEALLOCATE(m1)
       IF( ALLOCATED(n1) ) DEALLOCATE(n1)
       IF( ALLOCATED(rmomb) ) DEALLOCATE(rmomb)
       IF( ALLOCATED(zmomb) ) DEALLOCATE(zmomb)
       IF( ALLOCATED(mm) ) DEALLOCATE(mm)
       IF( ALLOCATED(dm1) ) DEALLOCATE(dm1)
       IF( ALLOCATED(faccon) ) DEALLOCATE(faccon)
       IF( ALLOCATED(xmpq) ) DEALLOCATE(xmpq)
       IF( ALLOCATED(result) ) DEALLOCATE(result)
       IF( ALLOCATED(rin3d) ) DEALLOCATE(rin3d)
       IF( ALLOCATED(zin3d) ) DEALLOCATE(zin3d)
       IF( ALLOCATED(angle) ) DEALLOCATE(angle)
!...end restore memory
       RETURN
       END
