C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C****************DIFFERENT-SUBROUTINE-USED-*****************************
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C		PRINT_SUBROUTINE
	
	SUBROUTINE PRINT4(F,I,LOUT)
	IMPLICIT NONE
	INTEGER MAXPOW,I0,IFLAG
	PARAMETER (MAXPOW = 10)
	INTEGER I, POWTEN, J, DIGIT(MAXPOW), K(MAXPOW), LOUT
	CHARACTER*2 F, A*10
	A='0123456789'
	IFLAG=0
	I0=10
	   DO J=1,I0
	   IF(IFLAG .EQ. 0)THEN
		IF((I/10**J).EQ.0)THEN
			POWTEN=J-1
			IFLAG=1
C	PAUSE
C	J=10
C	PRINT*,POWTEN,J,I0,I
		ENDIF
	      ENDIF
	   ENDDO
	   DO J=1,POWTEN+1
	DIGIT(J)=(I/10**(POWTEN+1-J))-I/10**(POWTEN-J+2)*10
	K(J)=DIGIT(J)+1
	   ENDDO
 60 	GOTO (10,20,30,40,50) POWTEN+1
 10	OPEN (UNIT=LOUT,FILE=F//A(1:1)//A(I+1:I+1),STATUS='UNKNOWN')
	GOTO 70
 20	OPEN(UNIT=LOUT,FILE=F//A(K(1):K(1))//A(K(2):K(2)),
     &	STATUS='UNKNOWN')
	GOTO 70
 30	OPEN (UNIT=LOUT,FILE=F//A(K(1):K(1))//A(K(2):K(2))//A(K(3):K(3))
     &	,STATUS='UNKNOWN')
	GOTO 70
 40	OPEN (UNIT=LOUT,FILE=F//A(K(1):K(1))//A(K(2):K(2))
     &	//A(K(3):K(3))//A(K(4):K(4)),STATUS='UNKNOWN')
	GOTO 70
 50	OPEN (UNIT=LOUT,FILE=F//A(K(1):K(1))//A(K(2):K(2))
     &	//A(K(3):K(3))//A(K(4):K(4))//A(K(5):K(5)),STATUS='UNKNOWN')
 70	RETURN
	END

