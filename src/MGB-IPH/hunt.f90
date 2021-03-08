	SUBROUTINE hunt(xx,x,jlo,yy,y,zz,z,n)
	!use nrtype
	implicit none
	integer, INTENT(INOUT) :: jlo
	real*8, INTENT(OUT) :: y,z
	real*8, INTENT(IN) :: x
	real*8, DIMENSION(n), INTENT(IN) :: xx, yy, zz
	integer :: n,inc,jhi,jm
	real*8 incjlo
	LOGICAL :: ascnd
	!n=size(xx)
	ascnd = (xx(n) >= xx(1))
	if (jlo <= 0 .or. jlo > n) then
		jlo=0
		jhi=n+1
	else
		inc=1
		if (x >= xx(jlo) .eqv. ascnd) then
			do
				jhi=jlo+inc
				if (jhi > n) then
					jhi=n+1
					exit
				else
					if (x < xx(jhi) .eqv. ascnd) exit
					jlo=jhi
					inc=inc+inc
				end if
			end do
		else
			jhi=jlo
			do
				jlo=jhi-inc
				if (jlo < 1) then
					jlo=0
					exit
				else
					if (x >= xx(jlo) .eqv. ascnd) exit
					jhi=jlo
					inc=inc+inc
				end if
			end do
		end if
	end if
	do
		if (jhi-jlo <= 1) then
			if (x == xx(n)) jlo=n-1
			if (x == xx(1)) jlo=1
			exit
		else
			jm=(jhi+jlo)/2
			if (x >= xx(jm) .eqv. ascnd) then
				jlo=jm
			else
				jhi=jm
			end if
		end if
	end do
	
	!Interpolação do vetor yy entre os valores jlo e jlo+1 (para Área)
	if(jlo<n)then
		incjlo = jlo + (1/(xx(jlo+1)-xx(jlo)))*(x-xx(jlo))
		y = yy(jlo) + (yy(jlo+1)-yy(jlo))*(incjlo-jlo)
	else
		y = yy(jlo) + (yy(jlo)-yy(jlo-5))/(5) * (jlo)
	endif
	
	!Interpolação do vetor yy entre os valores jlo e jlo+1 (para Nível)
	if(jlo<n)then
		incjlo = jlo + (1/(xx(jlo+1)-xx(jlo)))*(x-xx(jlo))
		z = zz(jlo) + (zz(jlo+1)-zz(jlo))*(incjlo-jlo)
	else
		z = zz(jlo) + (zz(jlo)-zz(jlo-5))/(5) * (jlo)
	endif
	
	end SUBROUTINE hunt
