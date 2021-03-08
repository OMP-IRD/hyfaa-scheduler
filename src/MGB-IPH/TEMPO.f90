	SUBROUTINE TEMPO(IBOBO,ISEED)
	!esta subrotina serve para estimar 
	!o tempo de processamento do programa
	implicit none
	integer IBOBO
	integer ITIME3,ITIME2,ITIME1 !   CONTADOR DE TEMPO
	integer ISEED
	
	TEMPOS: select case(IBOBO)
	case(0)   !   INICIA A CONTAGEM DE TEMPO
		itime1 = TIME()
		ISEED=ITIME1
		!print *, 'TEMPO INICIAL ', itime1
	case(1)	  !   TERMINA CONTAGEM DO TEMPO
		itime2 = TIME()
		ITIME3=ITIME2-ITIME1
		!print *, 'TEMPO TOTAL ', itime3,'SEG'
	end select TEMPOS
	RETURN
	end
