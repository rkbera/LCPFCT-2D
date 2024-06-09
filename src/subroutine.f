C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C**********(HLMHOLTZ SOLVER)TO CALCULATE THE VELOCITY*******************
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	SUBROUTINE HELM(SOURCE,OUTPUT,XLEFT,XRIGHT,YBOTOM,YTOP)


      INCLUDE 'param.i'

      INTEGER INX, INY, I, J
      INTEGER LDPHI, MBDCND, NBDCND, IERROR

      REAL*8 SOURCE(NPTX,NPTY), OUTPUT(NPTX,NPTY)
      REAL*8 BDA(NPTY+1), BDB(NPTY+1), BDC(NPTX+1), BDD(NPTX+1)
      REAL*8 XLEFT, XRIGHT, YBOTOM, YTOP
      REAL*8 COEFPHI, PERTRB
      REAL*8 W(NWORK)
      REAL*8 RHSH(NPTX+1,NPTY+1)

      INX = NPTX
      INY = NPTY
      LDPHI = NPTX+1
      COEFPHI = 0.0D0  
      MBDCND = 0
      NBDCND = 0
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	DO J = 1, NPTY
!	  DO I = 2, NPTX-1
!		RHSH(I,J) = PWXI(I,J) 
!	  ENDDO
!	ENDDO
!	DO I = 2, NPTX-1
! 		RHSH(I,NPTY+1)=RHSH(I,1)
!	ENDDO
!	DO J = 1, NPTY
!		RHSH(1,J)= VELXI(1,J)
!		RHSH(NPTX,J)= VELXI(NPTX,J)
!	ENDDO
!	RHSH(1,NPTY+1)=RHSH(1,1)
!	RHSH(NPTX,NPTY+1)=RHSH(NPTX,1)
!C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	DO I= 2,NPTX
	  DO J = 2,NPTY
		RHSH(I,J) = SOURCE(I,J) 
	  ENDDO
	ENDDO

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  DO J = 1,NPTY+1
		RHSH(1,J) = SOURCE(1,J) 
		RHSH(NPTX+1,J) = SOURCE(1,J) 
	  ENDDO
!C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	DO I=1,NPTX+1
		RHSH(I,1) = SOURCE(I,1) 
		RHSH(I,NPTY+1) = SOURCE(I,1) 
	ENDDO
	RHSH(1,NPTY+1)=RHSH(1,1)
	RHSH(NPTX+1,1)=RHSH(1,1)
	RHSH(NPTX+1,NPTY+1)=RHSH(1,NPTY+1)
!	RHSH(NPTX+1,NPTY+1)=RHSH(1,NPTY+1)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

         CALL HWSCRT(XLEFT,XRIGHT,INX,MBDCND,BDA,BDB,YBOTOM,YTOP,
     &        INY,NBDCND,BDC,BDD,COEFPHI,RHSH,LDPHI,PERTRB,IERROR,W)

      DO J = 1, NPTY
         DO I = 1, NPTX
            OUTPUT(I,J) = RHSH(I,J)
         ENDDO
      ENDDO

      RETURN
      END
C------------------------------------------------------




C SUBROUTINES FOR VELOCITY CALCULATION DERIVATIVES

C-------
	SUBROUTINE DERIVX(VALX,VALY,DX,DVALX,DVALY)
	INCLUDE 'param.i'
	INTEGER I,J
	INTEGER INX, INY
	REAL*8 VALX(NPTX,NPTY), VALY(NPTX,NPTY)
	REAL*8 DVALX(NPTX,NPTY), DVALY(NPTX,NPTY)
	REAL*8 DX
C-------
	INX = NPTX
	INY = NPTY
C-------
	   DO 11 J = 1,INY
	      DO 12 I = 2,INX-1
		DVALX(I,J) = (VALX(I+1,J) - VALX(I-1,J))/(2.0D0*DX)
		DVALY(I,J) = (VALY(I+1,J) - VALY(I-1,J))/(2.0D0*DX)
 12	      ENDDO
C-------
		DVALX(1,J)   = (VALX(2,J) - VALX(INX,J))/(2.0D0*DX)
		DVALX(INX,J) = (VALX(1,J) - VALX(INX-1,J))/(2.0D0*DX)
		DVALY(1,J)   = (VALY(2,J) - VALY(INX,J))/(2.0D0*DX)
		DVALY(INX,J) = (VALY(1,J) - VALY(INX-1,J))/(2.0D0*DX) 
C-------
!	   DVALX(1,J)   = DVALX(2,J)
!	   DVALX(INX,J) = DVALX(INX-1,J)
!	   DVALY(1,J)   = DVALY(2,J)
!	   DVALY(INX,J) = DVALY(INX-1,J) 
!C-------
 11	   ENDDO 
C-------
	RETURN
	END
C-------
C-----------------------------------------------------------------------
C	SUBROUTINE FOR DERIVATIVES IN Y DIRECTION WITH PBC AT CELL CENTRES
C-------
	SUBROUTINE DERIVY(VALX,VALY,DY,DVALX,DVALY)
C-------
	INCLUDE 'param.i'
	INTEGER I,J
	INTEGER INX, INY
	REAL*8 VALX(NPTX,NPTY), VALY(NPTX,NPTY)
	REAL*8 DVALX(NPTX,NPTY), DVALY(NPTX,NPTY)
	REAL*8 DY
C-------
	INX = NPTX
	INY = NPTY
C-------
	   DO 11 I = 1,INX
	      DO 12 J = 2,INY-1
		DVALX(I,J) = (VALX(I,J+1) - VALX(I,J-1))/(2.0D0*DY)
		DVALY(I,J) = (VALY(I,J+1) - VALY(I,J-1))/(2.0D0*DY)
 12	      ENDDO
C-------
		DVALX(I,1)   = (VALX(I,2) - VALX(I,INY))/(2.0D0*DY)
		DVALX(I,INY) = (VALX(I,1) - VALX(I,INY-1))/(2.0D0*DY)
		DVALY(I,1)   = (VALY(I,2) - VALY(I,INY))/(2.0D0*DY)
		DVALY(I,INY) = (VALY(I,1) - VALY(I,INY-1))/(2.0D0*DY) 
 11	   ENDDO 
C-------
	RETURN
	END
C-------














