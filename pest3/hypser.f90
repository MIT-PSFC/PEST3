	SUBROUTINE hypser(a,b,c,z,series,deriv)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	COMPLEX(DPC), INTENT(IN) :: a,b,c,z
	COMPLEX(DPC), INTENT(OUT) :: series,deriv
	INTEGER(I4B) :: n
	INTEGER(I4B), PARAMETER :: MAXIT=1000
	COMPLEX(DPC) :: aa,bb,cc,fac,temp
	deriv=cmplx(0.0_dp,0.0_dp,kind=dpc)
	fac=cmplx(1.0_dp,0.0_dp,kind=dpc)
	temp=fac
	aa=a
	bb=b
	cc=c
	do n=1,MAXIT
		fac=((aa*bb)/cc)*fac
		deriv=deriv+fac
		fac=fac*z/n
		series=temp+fac
		if (series == temp) RETURN
		temp=series
		aa=aa+1.0
		bb=bb+1.0
		cc=cc+1.0
	end do
	call nrerror('hypser: convergence failure')
	END SUBROUTINE hypser
