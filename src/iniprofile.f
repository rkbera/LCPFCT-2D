C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C***********************INITIAL_CONDITION******************************
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	SUBROUTINE INDATA(VELX,VELY,VELZ,DEN,VBLX,VBLY,VBLZ,DENB,
     &    ELECX,ELECY,ELECZ,BX,BY,BZ)
	INCLUDE 'param.i'
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	GRID VARIABLES 
	
	INTEGER NX, NY,NZ
	INTEGER NXP, NYP
	REAL*8 XCELL(NPTX), YCELL(NPTY)
	REAL*8 XLEFT,XRIGHT, YTOP,YBOTOM
	REAL*8 XINT(NXINT), YINT(NYINT)
	REAL*8 DX, DY
	REAL*8 LX, LY
	REAL*8 DEL, K
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	COMPONENTS OF VELOCITYPERTURBATION 
	
	REAL*8 VELX(NPTX,NPTY), VELY(NPTX,NPTY),VELZ(NPTX,NPTY)
	REAL*8 VBLX(NPTX,NPTY), VBLY(NPTX,NPTY),VBLZ(NPTX,NPTY)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	COMPONENTS OF DENSITY
	
	REAL*8 DEN(NPTX,NPTY),DENB(NPTX,NPTY),TOTALDEN(NPTX,NPTY)
	REAL*8 PHI(NPTX,NPTY),RHO(NPTX,NPTY)
	REAL*8 ELECY(NPTX,NPTY),PI
	REAL*8 ELEC(NPTX,NPTY),ELECX(NPTX,NPTY),ELECZ(NPTX,NPTY)
	REAL*8 BX(NPTX,NPTY), BY(NPTX,NPTY), BZ(NPTX,NPTY)
	REAL*8 DENX(NPTX),EXN(NPTX),EYN(NPTY),DENY(NPTY)
	REAL*8 JX(NPTX,NPTY), JY(NPTX,NPTY), DJXY(NPTX,NPTY)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	REAL*8 A,B,R(NPTX,NPTY),THETA
	 

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	DO LOOP INDICES
	INTEGER I, J
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	COMMON BLOCK
	COMMON / DATAP / NX,NY,NXP,NYP	
	COMMON / GRID / XINT,YINT,XCELL,YCELL,DX,DY,LX,LY,PI
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	NX = NPTX
	NY = NPTY
	
C	DENSITY+PERTURBATION 
	  DO  I = 1,NX
	  DO  J = 1,NY

	DEN(I,J)= 1.0d0 


c--------------------------------	
	
!	DENB(I,J)=7.0d0*DEXP(-((((XCELL(I)-5.0D0)**2)/2.4d0) 
!     &    +(YCELL(J)**2)/0.8D0)) 

!	if (XCELL(I) .ge. 1.0d0 .and. XCELL(I) .le. 2.0d0 ) then
!	if (YCELL(J) .ge. -2.0d0 .and. YCELL(J) .le. 2.0d0 ) then

!	  DENB(I,J)=0.5d0
	  DENB(I,J)=0.5d0*DEXP(-((((XCELL(I)-3.0D0)**2)/2.88d0) 
     &    +(YCELL(J)**2 )/0.32D0)) 

!	endif	
!	endif

	
	  TOTALDEN(I,J)=DEN(I,J)+DENB(I,J)
	
 	  ENDDO
 	  ENDDO

!C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!C	VELOCITY+PERTURBATION 

	DO 9 I = 1,NX
	   DO 10 J = 1,NY
		

		ELECX(I,J)=0.0D0
		ELECY(I,J)=0.0D0
		ELECZ(I,J)=0.0D0

		VELX(I,J)= 0.0D0	
		VELY(I,J)= 0.0D0
		VELZ(I,J)= 0.0D0

		VBLX(I,J)= 0.99999999D0	
		VBLY(I,J)= 0.0D0
		VBLZ(I,J)= 0.0D0

		BX(I,J)= 0.0D0
		BY(I,J)= 0.0D0

		BZ(I,J)= 0.0d0
 10	ENDDO
 9	ENDDO


c============================================================

!	
!	XLEFT  = XCELL(1)
!	XRIGHT = XCELL(NPTX)+DX
!	YBOTOM = YCELL(1)
!	YTOP   = YCELL(NPTY)+DY
!C-----------------------------------------
!	
!	DO  I = 1, NX
!	DO J = 1, NY
!				
!		JX(I,J)= -DENB(I,J)*VBLX(I,J)
!	 ENDDO
!	ENDDO

!	CALL DERIVY(JX,JX,DY,DJXY,DJXY)
!C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!C	TO CALCULATE VELOCITY USING HELMHOLTZ SOLVER
!C------	
!	CALL HELM(DJXY,BZ,XLEFT,XRIGHT,YBOTOM,YTOP)
!	

!C====================================================================
!	CALL PHI2D(TOTALDEN,ELECX,ELECY)
		
c--------------------------------------------------------------------------
	
	

	
c--------------------------------------------------
        RETURN
        END SUBROUTINE INDATA
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
