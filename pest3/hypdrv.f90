	SUBROUTINE hypdrv(s,ry,rdyds)
	USE nrtype
	USE hypgeo_info
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: s
	REAL(DP), DIMENSION(:), INTENT(IN) :: ry
	REAL(DP), DIMENSION(:), INTENT(OUT) :: rdyds
	COMPLEX(DPC), DIMENSION(2) :: y,dyds
	COMPLEX(DPC) :: z
	y=cmplx(ry(1:4:2),ry(2:4:2),kind=dpc)
	z=hypgeo_z0+s*hypgeo_dz
	dyds(1)=y(2)*hypgeo_dz
	dyds(2)=((hypgeo_aa*hypgeo_bb)*y(1)-(hypgeo_cc-&
		((hypgeo_aa+hypgeo_bb)+1.0_dp)*z)*y(2))*hypgeo_dz/(z*(1.0_dp-z))
	rdyds(1:4:2)=real(dyds)
	rdyds(2:4:2)=aimag(dyds)
	END SUBROUTINE hypdrv
