	SUBROUTINE pzextr(iest,xest,yest,yz,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: iest
	REAL(DP), INTENT(IN) :: xest
	REAL(DP), DIMENSION(:), INTENT(IN) :: yest
	REAL(DP), DIMENSION(:), INTENT(OUT) :: yz,dy
	INTEGER(I4B), PARAMETER :: IEST_MAX=16
	INTEGER(I4B) :: j,nv
	INTEGER(I4B), SAVE :: nvold=-1
	REAL(DP) :: delta,f1,f2
	REAL(DP), DIMENSION(size(yz)) :: d,tmp,q
	REAL(DP), DIMENSION(IEST_MAX), SAVE :: x
	REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE :: qcol
	nv=assert_eq(size(yz),size(yest),size(dy),'pzextr')
	if (iest > IEST_MAX) call &
		nrerror('pzextr: probable misuse, too much extrapolation')
	if (nv /= nvold) then
		if (allocated(qcol)) deallocate(qcol)
		allocate(qcol(nv,IEST_MAX))
		nvold=nv
	end if
	x(iest)=xest
	dy(:)=yest(:)
	yz(:)=yest(:)
	if (iest == 1) then
		qcol(:,1)=yest(:)
	else
		d(:)=yest(:)
		do j=1,iest-1
			delta=1.0_dp/(x(iest-j)-xest)
			f1=xest*delta
			f2=x(iest-j)*delta
			q(:)=qcol(:,j)
			qcol(:,j)=dy(:)
			tmp(:)=d(:)-q(:)
			dy(:)=f1*tmp(:)
			d(:)=f2*tmp(:)
			yz(:)=yz(:)+dy(:)
		end do
		qcol(:,iest)=dy(:)
	end if
	END SUBROUTINE pzextr
