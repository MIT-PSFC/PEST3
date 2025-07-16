	SUBROUTINE mmid(y,dydx,xs,htot,nstep,yout,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: nstep
	REAL(DP), INTENT(IN) :: xs,htot
	REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
	REAL(DP), DIMENSION(:), INTENT(OUT) :: yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: n,ndum
	REAL(DP) :: h,h2,x
	REAL(DP), DIMENSION(size(y)) :: ym,yn
	ndum=assert_eq(size(y),size(dydx),size(yout),'mmid')
	h=htot/nstep
	ym=y
	yn=y+h*dydx
	x=xs+h
	call derivs(x,yn,yout)
	h2=2.0_dp*h
	do n=2,nstep
		call swap(ym,yn)
		yn=yn+h2*yout
		x=x+h
		call derivs(x,yn,yout)
	end do
	yout=0.5_dp*(ym+yn+h*yout)
	END SUBROUTINE mmid
