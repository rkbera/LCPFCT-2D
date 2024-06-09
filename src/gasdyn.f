C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C***************INCOMPRESIBLE-GHD-EVOLUTION*******************
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	SUBROUTINE GASDYN
C------
	SUBROUTINE GASD1(DT,VELX,VELY,VELZ,DEN,VBLX,VBLY,VBLZ,DENB,
     &   EX,EY,EZ,BX,BY,BZ,TSTEP) 


C------
C	DECLARATION BLOCK
	
	INCLUDE 'param.i'
C------
C	DO LOOP INDICES
	
	INTEGER I, J 
	INTEGER IT, TSTEP
C------
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	GRID VARIABLES
C------
	INTEGER I1, INX, INY,INZ
	INTEGER INXP, INYP ,INZP
	REAL*8  XINT(NXINT), YINT(NYINT),ZINT(NYINT)
	REAL*8  XCELL(NPTX), YCELL(NPTY),ZCELL(NPTY)
	REAL*8  DX, DY
	REAL*8  LX, LY
	REAL*8  XLEFT,XRIGHT,YTOP,YBOTOM

C------
C==============================================================
C	********************************************
C 	*                                         *
C 	*					  *	
C 	*	BACKGROUND PLASMA PARAMETERS      *
C 	*                                         *
C	********************************************
C=============================================================
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	COMPONENTS OF PLASMA VELOCITY
C------

	REAL*8 PEX(NPTX,NPTY), PEY(NPTX,NPTY),PEZ(NPTX,NPTY)
	REAL*8 PEXI(NPTX,NPTY),PEYI(NPTX,NPTY), PEZI(NPTX,NPTY)
C------
	REAL*8 PEX0(NPT),PEY0(NPT),PEZ0(NPT)
	
	REAL*8 PEXN(NPT),PEYN(NPT), PEZN(NPT)
	
	REAL*8 PE(NPTX,NPTY)
	REAL*8 PEI(NPTX,NPTY)
C------
	REAL*8 PE0(NPT)
	
	REAL*8 PEN(NPT)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	REAL*8 VELX(NPTX,NPTY), VELY(NPTX,NPTY),VELZ(NPTX,NPTY)
	REAL*8 VELXI(NPTX,NPTY),VELYI(NPTX,NPTY), VELZI(NPTX,NPTY)
C------
	REAL*8 VELX0(NPT),VELY0(NPT),VELZ0(NPT)
	REAL*8 VINT(NPT) 
	REAL*8 VELXN(NPT),VELYN(NPT), VELZN(NPT)
C---------
C	COMPONENTS OF PLASMA DENSITY 
	
	REAL*8 DENX(NPTX,NPTY),DENY(NPTX,NPTY),DENZ(NPTX,NPTY)
	REAL*8 DENXI(NPTX,NPTY),DENYI(NPTX,NPTY), DENZI(NPTX,NPTY)

	REAL*8 DEN(NPTX,NPTY), DEN0(NPTX,NPTY)
	REAL*8 DENI(NPTX,NPTY)
	
	REAL*8 DENX0(NPT),DENXN(NPT)
	REAL*8 DENY0(NPT),DENYN(NPT)
	REAL*8 DENZ0(NPT),DENZN(NPT)
c---------

C       COMPONENTS FOR ELECTRIC FIELD

	REAL*8 EX(NPTX,NPTY),EY(NPTX,NPTY), EZ(NPTX,NPTY)
	REAL*8 EXI(NPTX,NPTY),EYI(NPTX,NPTY), EZI(NPTX,NPTY)

	REAL*8 EX0(NPT),EY0(NPT),EZ0(NPT)
	REAL*8 EXN(NPT),EYN(NPT),EZN(NPT)

C-----------------------------------------	
	REAL*8 PHIX(NPTX,NPTY),PHIY(NPTX,NPTY)
	REAL*8 PHIXI(NPTX,NPTY),PHIYI(NPTX,NPTY)

	REAL*8 PHI(NPTX,NPTY), PHI0(NPTX,NPTY)
	REAL*8 PHII(NPTX,NPTY)
	
	REAL*8 PHIX0(NPT),PHIXN(NPT)
	REAL*8 PHIY0(NPT),PHIYN(NPT)
	
c---------
C-------------------------------------------------------

C	COMPONENTS FOR CURRENT DENSITY

	

	REAL*8 JX(NPTX,NPTY),JY(NPTX,NPTY), JZ(NPTX,NPTY)
	REAL*8 JXI(NPTX,NPTY),JYI(NPTX,NPTY), JZI(NPTX,NPTY)

	REAL*8 JX0(NPT),JY0(NPT),JZ0(NPT)
	REAL*8 JXN(NPT),JYN(NPT),JZN(NPT)


C       COMPONENTS FOR MAGNETIC FIELD

	REAL*8 BX(NPTX,NPTY),BY(NPTX,NPTY), BZ(NPTX,NPTY)
	REAL*8 BXI(NPTX,NPTY),BYI(NPTX,NPTY), BZI(NPTX,NPTY)

	REAL*8 BX0(NPT),BY0(NPT),BZ0(NPT)
	REAL*8 BXN(NPT),BYN(NPT),BZN(NPT)

C---------
C	COMPONENTS FOR RELATIVISTIC FACTOR
	
	REAL*8 GAMAX(NPTX,NPTY),GAMAY(NPTX,NPTY),GAMAZ(NPTX,NPTY)
	REAL*8 GAMAXI(NPTX,NPTY),GAMAYI(NPTX,NPTY),GAMAZI(NPTX,NPTY)

	REAL*8 GAMA(NPTX,NPTY), GAMA0(NPTX,NPTY)
	REAL*8 GAMAI(NPTX,NPTY),GAMAN(NPT)
	
	REAL*8 GAMAX0(NPT),GAMAXN(NPT)
	REAL*8 GAMAY0(NPT),GAMAYN(NPT)	
	REAL*8 GAMAZ0(NPT),GAMAZN(NPT)




C	SOURCE VARIABLES FOR PLASMA X DIRECTION
	
	REAL*8 CXX1(NPTX), CXX2(NPTX),CXX3(NPTX),CXX4(NPTX)
	REAL*8 CXX5(NPTX), CXX6(NPTX),CXX7(NPTX)
	REAL*8 CXX8(NPTX), CXX9(NPTX),CXX10(NPTX)
	REAL*8 CXX11(NPTX), CXX12(NPTX),CXX13(NPTX),CXX14(NPTX)

	REAL*8 DXX1(NPTX),DXX2(NPTX),DXX3(NPTX),DXX4(NPTX)
	REAL*8 DXX5(NPTX),DXX6(NPTX),DXX7(NPTX)
	REAL*8 DXX11(NPTX),DXX12(NPTX),DXX13(NPTX),DXX14(NPTX)
	REAL*8 DXX8(NPTX),DXX9(NPTX),DXX10(NPTX)

	REAL*8 DXX1I1,DXX2I1,DXX3I1,DXX4I1,DXX5I1,DXX6I1,DXX7I1
	REAL*8 DXX8I1,DXX9I1,DXX10I1,DXX11I1,DXX12I1,DXX13I1,DXX14I1

	REAL*8 DXX1INP,DXX2INP,DXX3INP,DXX4INP,DXX5INP,DXX6INP
	REAL*8 DXX11INP,DXX12INP,DXX13INP,DXX14INP
	REAL*8 DXX7INP, DXX8INP,DXX9INP,DXX10INP

c----------------------------------------------------------------
C	SOURCE VARIABLES  Y DIRECTION
	
	REAL*8 CYY1(NPTY), CYY2(NPTY),CYY3(NPTY),CYY4(NPTY)
	REAL*8 CYY5(NPTY), CYY6(NPTY),CYY7(NPTY)
	REAL*8 CYY8(NPTY), CYY9(NPTY),CYY10(NPTY)
	REAL*8 CYY11(NPTY), CYY12(NPTY),CYY13(NPTY),CYY14(NPTY)
	REAL*8 DYY1(NPTY), DYY2(NPTY),DYY3(NPTY), DYY4(NPTY)
	REAL*8 DYY5(NPTY), DYY6(NPTY),DYY7(NPTY)
	REAL*8 DYY8(NPTY), DYY9(NPTY),DYY10(NPTY)
	REAL*8 DYY11(NPTY), DYY12(NPTY),DYY13(NPTY),DYY14(NPTY)
	REAL*8 DYY1I1,DYY2I1,DYY3I1,DYY4I1,DYY5I1,DYY6I1,DYY7I1
	REAL*8 DYY8I1,DYY9I1,DYY10I1,DYY11I1,DYY12I1,DYY13I1,DYY14I1

	REAL*8 DYY1INP,DYY2INP,DYY3INP,DYY4INP,DYY5INP,DYY6INP

	REAL*8 DYY7INP,DYY8INP,DYY9INP,DYY10INP,DYY11INP
	REAL*8 DYY12INP,DYY13INP,DYY14INP
C	


C     COMPONENTS OF BALANCING FORCE

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	VARIABLES OF THE BOUNDARY CONDITIONS OF PLASMA DENSITY
C------
	REAL*8 SRHO1X, SRHONX, VRHO1X, VRHONX
	REAL*8 SRHO1Y, SRHONY, VRHO1Y, VRHONY
	REAL*8 SRHO1Z, SRHONZ, VRHO1Z, VRHONZ
C==============================================================
C	********************************************
C 	*                                          *
C 	*					   *	
C 	*	BEAM PARAMETERS                    *
C 	*                                          *
C	********************************************
C=============================================================
	REAL*8 VBLX(NPTX,NPTY), VBLY(NPTX,NPTY), VBLZ(NPTX,NPTY)
	REAL*8 VBLXI(NPTX,NPTY),VBLYI(NPTX,NPTY), VBLZI(NPTX,NPTY)
C------------------------------------------------------------
	REAL*8 VBLX0(NPT),VBLY0(NPT),VBLZ0(NPT)
	REAL*8 VBINT(NPT) 
	REAL*8 VBLXN(NPT),VBLYN(NPT),VBLZN(NPT)
C------------------------------------------------------------
C	COMPONENTS OF BEAM DENSITY 
	
	REAL*8 DENBX(NPTX,NPTY),DENBY(NPTX,NPTY),DENBZ(NPTX,NPTY)
	REAL*8 DENBXI(NPTX,NPTY),DENBYI(NPTX,NPTY), DENBZI(NPTX,NPTY)

	REAL*8 DENB(NPTX,NPTY), DENB0(NPTX,NPTY)
	REAL*8 DENBI(NPTX,NPTY)
	
	REAL*8 DENBX0(NPT),DENBXN(NPT)
	REAL*8 DENBY0(NPT),DENBYN(NPT)
	REAL*8 DENBZ0(NPT),DENBZN(NPT)
C	
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	COMPONENTS OF BEAM VELOCITY
C------

	REAL*8 PBX(NPTX,NPTY), PBY(NPTX,NPTY),PBZ(NPTX,NPTY)
	REAL*8 PBXI(NPTX,NPTY),PBYI(NPTX,NPTY), PBZI(NPTX,NPTY)
C------
	REAL*8 PBX0(NPT),PBY0(NPT),PBZ0(NPT)
	
	REAL*8 PBXN(NPT),PBYN(NPT), PBZN(NPT)
	
	REAL*8 PB(NPTX,NPTY)
	REAL*8 PBI(NPTX,NPTY)
C------
	REAL*8 PB0(NPT)
	
	REAL*8 PBN(NPT)

C-------------------------
C	COMPONENTS FOR RELATIVISTIC FACTOR
	
	REAL*8 GAMABX(NPTX,NPTY),GAMABY(NPTX,NPTY),GAMABZ(NPTX,NPTY)
	REAL*8 GAMABXI(NPTX,NPTY),GAMABYI(NPTX,NPTY),GAMABZI(NPTX,NPTY)

	REAL*8 GAMAB(NPTX,NPTY), GAMAB0(NPTX,NPTY)
	REAL*8 GAMABI(NPTX,NPTY),GAMABN(NPT)
C--------------------------------------------------------------
C	FULL AND HALF TIME STEPS
	
	REAL*8 DT, DTSUB
	
C	LOGICAL VARIABLE TO DETERMINE THE TYPE OF BOUNDARY CONDITIONS 
	
	LOGICAL PBC
	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	COMMON-GLOBLE-BLOK
C------
	COMMON /GRID/ XINT,YINT,XCELL,YCELL,DX,DY,LX,LY
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	INITIALIZING VARIABLES
C------
	I1 = 1
	INX = NPTX
	INY = NPTY
C------
	INXP = INX + 1
	INYP = INY + 1	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	INITIALIZATION OF INTERMEDIATE 2-D  ARRAYS
C------
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	MOMENTUM CALCULATION

	DO I=1,INX
	DO J=1,INY
	GAMA(I,J)=1.0D0/DSQRT(1-((VELX(I,J))**2+(VELY(I,J))**2+
     &    (VELZ(I,J))**2))

	PEX(I,J)= GAMA(I,J)*VELX(I,J)
	
	PEY(I,J)= GAMA(I,J)*VELY(I,J)

	PEZ(I,J)= GAMA(I,J)*VELZ(I,J)		

	ENDDO
	ENDDO




C------------------------------------------------------------------

	DO 2 I = 1, INX
	   DO 3 J = 1, INY

		DENI(I,J)= DEN(I,J)

		VELXI(I,J) = VELX(I,J)
		VELYI(I,J) = VELY(I,J)
		VELZI(I,J) = VELZ(I,J)

		PEXI(I,J) = PEX(I,J)
		PEYI(I,J) = PEY(I,J)
		PEZI(I,J) = PEZ(I,J)



		DENBI(I,J)= DENB(I,J)

		VBLXI(I,J) = VBLX(I,J)
		VBLYI(I,J) = VBLY(I,J)
		VBLZI(I,J) = VBLZ(I,J)


		
		EXI(I,J)= EX(I,J)
		EYI(I,J)= EY(I,J)
		EZI(I,J)= EZ(I,J)
		
		BXI(I,J)= BX(I,J)
		BYI(I,J)= BY(I,J)
		BZI(I,J)= BZ(I,J)
		
3	   CONTINUE
2	CONTINUE
C------
c
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C***************INTEGRATION IN X-DIRECTION STARTS HERE *****************
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	PBC = .FALSE.
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	PARAMETERS FOR THE BOUNDARY CONDITIONS OF VORTICITY
C------
	SRHO1X = -1.0D0
	SRHONX = -1.0D0
	VRHO1X =  0.0D0
	VRHONX =  0.0D0
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	SETTING GRID FACTORS FOR USE IN CNVFCT
	
	CALL RESIDIFF( 1.0D0 )
	CALL MAKEGRID( XINT, XINT, I1, INXP, 1 )
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C------
	DO 90 IT = 1,2
		DTSUB = 0.5D0*DFLOAT(IT)*DT
C------
	
	        DO 80 J = 1, INY
		DO 70 I = 1, INX 

		DENX0(I) = DEN(I,J)  
  
		DENBX0(I) = DENB(I,J)

		VELX0(I)=VELX(I,J)
		VELY0(I)=VELY(I,J)
		VELZ0(I)=VELZ(I,J)	
		
		PEX0(I)=PEX(I,J)
		PEY0(I)=PEY(I,J)
		PEZ0(I)=PEZ(I,J)
		 
			

		EX0(I)= EX(I,J)
		EY0(I)= EY(I,J)
		EZ0(I)= EZ(I,J)

		BX0(I)= BX(I,J)
		BY0(I)= BY(I,J)
		BZ0(I)= BZ(I,J)
		
		
C-----------------------------------------------------------------------
C	HERE WE USE THE SOURCE TERM FOR LCPFCT IN X
C-----------------------------------------------------------------------


C------VELOCITY FOR PLASMA --------------------------------

	CXX1(I) = 1.0D0
	DXX1(I) = -EXI(I,J)-VELYI(I,J)*BZI(I,J)+VELZI(I,J)*BYI(I,J)


	CXX2(I) = 1.0D0
	DXX2(I) = -EYI(I,J)-VELZI(I,J)*BXI(I,J)+VELXI(I,J)*BZI(I,J)

	CXX3(I) = 1.0D0
	DXX3(I) = -EZI(I,J)-VELXI(I,J)*BYI(I,J)+VELYI(I,J)*BXI(I,J)


C-------ELECTRIC FIELD

	CXX4(I) = 1.0D0
	DXX4(I) = DENI(I,J)*VELXI(I,J)+ DENBI(I,J)*VBLXI(I,J)


	CXX5(I) = 1.0D0
	DXX5(I) = DENI(I,J)*VELYI(I,J)+ DENBI(I,J)*VBLYI(I,J)

	CXX6(I) = 1.0D0
	DXX6(I) = -BZI(I,J)

	CXX7(I) = 1.0D0
	DXX7(I) = DENI(I,J)*VELZI(I,J)+ DENBI(I,J)*VBLZI(I,J)

	
	CXX8(I) = 1.0D0
	DXX8(I) = BYI(I,J)

C----- MGNETIC FIELD

	CXX9(I) = 1.0D0
	DXX9(I) = 0.0D0

	CXX10(I) = 1.0D0
	DXX10(I) = EZI(I,J)

		
	CXX11(I) = 1.0D0
	DXX11(I) = -EYI(I,J)



C-----------------------------------------------------------------------
 70		CONTINUE
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C------ DENSITY

	DXX1I1= (DXX1(1)+DXX1(INX))/(2.0D0)         
	DXX1INP = DXX1(INX)	

	DXX2I1= (DXX2(1)+DXX2(INX))/(2.0D0) 
	DXX2INP = DXX2(INX)	
	
	DXX3I1= (DXX3(1)+DXX3(INX))/(2.0D0) 
	DXX3INP = DXX3(INX)
	
	DXX4I1= (DXX4(1)+DXX4(INX))/(2.0D0) 
	DXX4INP = DXX4(INX)	

	DXX5I1= (DXX5(1)+DXX5(INX))/(2.0D0) 
	DXX5INP = DXX5(INX)	
	
	DXX6I1= (DXX6(1)+DXX6(INX))/(2.0D0) 
	DXX6INP = DXX6(INX)

	DXX7I1= (DXX7(1)+DXX7(INX))/(2.0D0) 
	DXX7INP = DXX7(INX)

	DXX8I1= (DXX8(1)+DXX8(INX))/(2.0D0) 
	DXX8INP = DXX8(INX)

	DXX9I1= (DXX9(1)+DXX9(INX))/(2.0D0) 
	DXX9INP = DXX9(INX)

	DXX10I1= (DXX10(1)+DXX10(INX))/(2.0D0) 
	DXX10INP = DXX10(INX)	

	DXX11I1= (DXX11(1)+DXX11(INX))/(2.0D0) 
	DXX11INP = DXX11(INX)



C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CC	SETTING X-INTERFACE VELOCITIES
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C	PLASMA EVOLUTION
	
		DO 30 I = I1 + 1, INX
		    VINT(I) = (VELXI(I-1,J) + VELXI(I,J))/2.0D0
		  
 30		CONTINUE
	
C	SETTING X-BOUNDRIES INTERFACE VELOCITIES
	
		VINT(I1) = (VELXI(1,J)+VELXI(INX,J))/2.0D0
		VINT(INXP) = VINT(INX)

		
C	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	CALLING SUBROUTINES VELOCITY,SOURCES AND CNVFCT
	
		CALL VELOCITY(VINT, I1, INXP, DTSUB)
		

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C		PLASMA ELECTRON EVOLUTION


		CALL LCPFCT(DENX0,DENXN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C		PLASMA VELOCITY ALONG X DIRECTION EVOLUTION
		
	        CALL SOURCES(I1,INX,DTSUB,3,CXX1,DXX1,DXX1I1,DXX1INP)
	

		CALL CNVFCT(PEX0,PEXN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)


C--------------------------------------------------------------

		CALL SOURCES(I1,INX,DTSUB,3,CXX2,DXX2,DXX2I1,DXX2INP)


		CALL CNVFCT(PEY0,PEYN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)
	


		CALL SOURCES(I1,INX,DTSUB,3,CXX3,DXX3,DXX3I1,DXX3INP)

		CALL CNVFCT(PEZ0,PEZN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)


C====================================================================
C	FIELD EVOLUTION
C------------------------------------------------------------------
	
		DO  I = I1 + 1, INX
		    VINT(I) = 0.0D0
		ENDDO
	
C	SETTING X-BOUNDRIES INTERFACE VELOCITIES
	
		VINT(I1) = 0.0D0
		VINT(INXP) = VINT(INX)

		
C	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	CALLING SUBROUTINES VELOCITY,SOURCES AND CNVFCT
	
		CALL VELOCITY(VINT, I1, INXP, DTSUB)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C		ELECTRIC FIELD COMPONENT

		CALL SOURCES(I1,INX,DTSUB,3,CXX4,DXX4,DXX4I1,DXX4INP)


		CALL CNVFCT(EX0,EXN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)


C------------------------------------------------------------------


		CALL SOURCES(I1,INX,DTSUB,3,CXX5,DXX5,DXX5I1,DXX5INP)

  		CALL SOURCES(I1,INX,DTSUB,2,CXX6,DXX6,DXX6I1,DXX6INP)

		
		CALL CNVFCT(EY0,EYN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)

C---------------------
		CALL SOURCES(I1,INX,DTSUB,3,CXX7,DXX7,DXX7I1,DXX7INP)

		CALL SOURCES(I1,INX,DTSUB,2,CXX8,DXX8,DXX8I1,DXX8INP)


		CALL CNVFCT(EZ0,EZN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	MAGNETIC FIELD COMPONENTS

		CALL SOURCES(I1,INX,DTSUB,3,CXX9,DXX9,DXX9I1,DXX9INP)


		CALL CNVFCT(BX0,BXN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)



		CALL SOURCES(I1,INX,DTSUB,2,CXX10,DXX10,DXX10I1,DXX10INP)
		
		CALL CNVFCT(BY0,BYN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)


		CALL SOURCES(I1,INX,DTSUB,2,CXX11,DXX11,DXX11I1,DXX11INP)

		
		CALL CNVFCT(BZ0,BZN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C 		BEAM DENSITY EVOLUTION
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		DO 31 I = I1 + 1, INX

		    VBINT(I) = (VBLXI(I-1,J) + VBLXI(I,J))/2.0D0
		  
 31		CONTINUE

		VBINT(I1) = (VBLXI(1,J)+VBLXI(INX,J))/2.0D0

		VBINT(INXP) = VBINT(INX)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C		BEAM EVOLUTION EQUATION 

		CALL VELOCITY( VBINT, I1, INXP, DTSUB)


		CALL LCPFCT(DENBX0,DENBXN,I1,INX,SRHO1X,VRHO1X,SRHONX,
     &		VRHONX,PBC)
C
C------------------------------------------------------------------------
C	RESULT AFTER EACH HALF TIME STEP

	   DO 40 I = 1, INX
		DENI(I,J) = DENXN(I)
		DENBI(I,J) = DENBXN(I)

		PEXI(I,J) = PEXN(I)
		PEYI(I,J) = PEYN(I)
		PEZI(I,J) = PEZN(I)


		
		EXI(I,J)= EXN(I)
		EYI(I,J)= EYN(I)
		EZI(I,J) = EZN(I)

		BXI(I,J)= BXN(I)
		BYI(I,J)= BYN(I)
		BZI(I,J)= BZN(I)
		
40	   CONTINUE

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	VELOCITY FORM MOMENTUM
	
	
	DO I=1,INX
	
	GAMAI(I,J)=DSQRT(1.0D0+((PEXI(I,J))**2+(PEYI(I,J))**2+
     &    (PEZI(I,J))**2))

	VELXI(I,J)= PEXI(I,J)/GAMAI(I,J)
	
	VELYI(I,J)= PEYI(I,J)/GAMAI(I,J)

	VELZI(I,J)= PEZI(I,J)/GAMAI(I,J)	

	
	ENDDO

C	VELOCITY OD SECOND BEAM FROM MOMENTUM
	
	





C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
80	CONTINUE          !CONTINUATION OF LOOP IN Y-DIRECTION 
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
90	CONTINUE                 ! END OF HALF-FULL TIME STEP	
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	FINAL RESULT AFTER COMPLITION OF TIME STEP
	   DO 1000 J = 1, INY
	      DO 1001 I = 1, INX

		DEN(I,J) = DENI(I,J)
		DENB(I,J) = DENBI(I,J)

		VELX(I,J) = VELXI(I,J)
		VELY(I,J) = VELYI(I,J)
		VELZ(I,J) = VELZI(I,J)


		PEX(I,J) = PEXI(I,J)
		PEY(I,J) = PEYI(I,J)
		PEZ(I,J) = PEZI(I,J)


		EX(I,J) = EXI(I,J)
		EY(I,J) = EYI(I,J)
		EZ(I,J) = EZI(I,J)

		BX(I,J) = BXI(I,J)
		BY(I,J) = BYI(I,J)
		BZ(I,J) = BZI(I,J)
		
 1001	      CONTINUE
 1000	   CONTINUE
!C
C
C
C********************************************************************	
C*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
C*	INTEGRATION ALONG Y-DIRECTION STARTS HERE
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	SETTING LOGICAL VARIABLE FOR PERIODIC BOUNDARY CONDITIONS IN Y

C*********************************************************************
C
C
C	
	PBC = .FALSE.
	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	PARAMETERS FOR THE BOUNDARY CONDITIONS OF DEN
C------
	SRHO1Y = -1.0D0
	SRHONY = -1.0D0
	VRHO1Y =  0.0D0
	VRHONY =  0.0D0

	
C------
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	SETTING GRID FACTORS FOR USE IN CNVFCT
C------
	CALL RESIDIFF( 1.0D0 )
	CALL MAKEGRID( YINT, YINT, I1, INYP, 1 )
C------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	BEGINING OF HALF-FULL TIME STEP INTEGRATION IN Y-DIRECTION
C------
	DO 900 IT = 1,2
		DTSUB = 0.5D0*DFLOAT(IT)*DT

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C------
	   DO 800 I = 1, INX
	      DO 700 J = 1, INY 
		DENY0(J)= DEN(I,J)        ! HERE WE USE DEN NOT DENI
		 
		VELX0(J)= VELX(I,J)
	        VELY0(J)= VELY(I,J)
		VELZ0(J)= VELZ(I,J)

		PEX0(J)= PEX(I,J)
	        PEY0(J)= PEY(I,J)
		PEZ0(J)= PEZ(I,J)

		DENBY0(J)= DENB(I,J)        ! HERE WE USE DEN NOT DENI
		 
		
		EX0(J) = EX(I,J)
		EY0(J)=EY(I,J)	
		EZ0(J)=EZ(I,J)

		BX0(J) = BX(I,J)
		BY0(J)=BY(I,J)	
		BZ0(J)=BZ(I,J)	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C-------VELOCITY SOURCE TERMS
	CYY1(J) = 1.0D0
	DYY1(J) =  0.0D0 

	
	CYY2(J) = 1.0D0
	DYY2(J) =  0.0D0 

	CYY3(J) = 1.0D0
	DYY3(J) =0.0D0 

C-------ELECTRIC FIELD COMPONENTS

	CYY4(J) = 1.0D0
	DYY4(J) =  0.0D0 

	CYY5(J) = 1.0D0
	DYY5(J) = BZI(I,J)

	CYY6(J) = 1.0D0
	DYY6(J) =  0.0D0 

	CYY7(J) = 1.0D0
	DYY7(J) = 0.0D0 
	
	CYY8(J) = 1.0D0
	DYY8(J) = -BXI(I,J)

C---- MGNETIC FIELD COMPONENTS

	CYY9(J) = 1.0D0
	DYY9(J) = -EZI(I,J)

	CYY10(J) = 1.0D0
	DYY10(J) = 0.0D0

	CYY11(J) = 1.0D0
	DYY11(J) = EXI(I,J)

	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 700		CONTINUE
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C------
C------DENSITY

	DYY1I1= (DYY1(1)+DYY1(INX))/(2.0D0)
	DYY1INP = DYY1(INY)
	

	DYY2I1= (DYY2(1)+DYY2(INX))/(2.0D0)
	DYY2INP = DYY2(INY)


	DYY3I1= (DYY3(1)+DYY3(INX))/(2.0D0)
	DYY3INP = DYY3(INY)

	DYY4I1= (DYY4(1)+DYY4(INX))/(2.0D0)
	DYY4INP = DYY4(INY)
	

	DYY5I1= (DYY5(1)+DYY5(INX))/(2.0D0)
	DYY5INP = DYY5(INY)


	DYY6I1= (DYY6(1)+DYY6(INX))/(2.0D0)
	DYY6INP =DYY6(INY)
	
	DYY7I1= (DYY7(1)+DYY7(INX))/(2.0D0)
	DYY7INP = DYY7(INY)

	DYY8I1= (DYY8(1)+DYY8(INX))/(2.0D0)
	DYY8INP = DYY8(INY)

	DYY9I1= (DYY9(1)+DYY9(INX))/(2.0D0)
	DYY9INP = DYY9(INX)

	DYY10I1= (DYY10(1)+DYY10(INX))/(2.0D0)
	DYY10INP = DYY10(INY)

	DYY11I1= (DYY11(1)+DYY11(INX))/(2.0D0)
	DYY11INP = DYY11(INY)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
C	SETTING Y-INTERFACE VELOCITIES
	
	   DO 300 J = I1 + 1, INY

		VINT(J) = (VELYI(I,J-1) + VELYI(I,J))/2.0D0
		
 300	   CONTINUE
	
C	SETTING Y-BOUNDRIES INTERFACE VELOCITIES
	
		VINT(I1) = (VELYI(I,1) + VELYI(I,INY))/2.0D0

		VINT(INYP) = VINT(I1)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	CALLING SUBROUTINES VELOCITY,SOURCES AND CNVFCT
C------
	CALL VELOCITY( VINT, I1, INYP, DTSUB )
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C    	 PLASMA---------------------
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	

C
	CALL LCPFCT(DENY0,DENYN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
C	PLASMA VELOCITY ALONG Y DIRECTION 

	CALL SOURCES(I1,INY,DTSUB,3,CYY1,DYY1,DYY1I1,DYY1INP)

	CALL CNVFCT(PEX0,PEXN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)
	
	        
	CALL SOURCES(I1,INY,DTSUB,3,CYY2,DYY2,DYY2I1,DYY2INP)
	
	CALL CNVFCT(PEY0,PEYN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)

	CALL SOURCES(I1,INY,DTSUB,3,CYY3,DYY3,DYY3I1,DYY3INP)
	
	CALL CNVFCT(PEZ0,PEZN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)
C=======================================================================
C	FIELD EVOLUTION
C------------------------------------------
	DO J = I1 + 1, INY

		VINT(J) = 0.0D0
		
 	ENDDO
	
C	SETTING Y-BOUNDRIES INTERFACE VELOCITIES
	
		VINT(I1) = 0.0D0

		VINT(INYP) = VINT(I1)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	CALLING SUBROUTINES VELOCITY,SOURCES AND CNVFCT
C------
	CALL VELOCITY( VINT, I1, INYP, DTSUB )
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	ELECTRIC FIELD COMPONENTS
	
	     
	CALL SOURCES(I1,INY,DTSUB,3,CYY4,DYY4,DYY4I1,DYY4INP)

	CALL SOURCES(I1,INY,DTSUB,2,CYY5,DYY5,DYY5I1,DYY5INP)
	

	CALL CNVFCT(EX0,EXN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)
   

	CALL SOURCES(I1,INY,DTSUB,3,CYY6,DYY6,DYY6I1,DYY6INP)
	
	CALL CNVFCT(EY0,EYN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)

	
	CALL SOURCES(I1,INY,DTSUB,3,CYY7,DYY7,DYY7I1,DYY7INP)

	CALL SOURCES(I1,INY,DTSUB,2,CYY8,DYY8,DYY8I1,DYY8INP)
	
	CALL CNVFCT(EZ0,EZN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	MAGNETIC FIELD COMPONENTS
	
	     

	CALL SOURCES(I1,INY,DTSUB,2,CYY9,DYY9,DYY9I1,DYY9INP)
	
	
	CALL CNVFCT(BX0,BXN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)
   

	
	CALL SOURCES(I1,INY,DTSUB,3,CYY10,DYY10,DYY10I1,DYY10INP)

	CALL CNVFCT(BY0,BYN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)

	
	CALL SOURCES(I1,INY,DTSUB,2,CYY11,DYY11,DYY11I1,DYY11INP)
	
	CALL CNVFCT(BZ0,BZN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C 	BEAM EVOLUTION IN Y

	DO 3001 J = I1 + 1, INY

		VBINT(J) = (VBLYI(I,J-1) + VBLYI(I,J))/2.0D0
3001 	continue
	
C	SETTING Y-BOUNDRIES INTERFACE VELOCITIES
	
		VBINT(I1) = (VBLYI(I,1) + VBLYI(I,INY))/2.0D0

		VBINT(INYP) = VBINT(INY)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	CALLING SUBROUTINES VELOCITY,SOURCES AND CNVFCT
C------
	CALL VELOCITY( VBINT, I1, INYP, DTSUB )
C--------------------------------------------------------------------

	CALL LCPFCT(DENBY0,DENBYN,I1,INY,SRHO1Y,VRHO1Y,SRHONY,
     &		VRHONY,PBC)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


C-----------------------------------------------------------

	   DO 400 J = 1, INY

		DENI(I,J)  = DENYN(J)
		DENBI(I,J)  = DENBYN(J)

		PEXI(I,J)  = PEXN(J)
		PEYI(I,J)  = PEYN(J)
		PEZI(I,J)  = PEZN(J)


		EXI(I,J) = EXN(J)
		EYI(I,J)= EYN(J)
		EZI(I,J)= EZN(J)
		
		BXI(I,J) = BXN(J)
		BYI(I,J)= BYN(J)
		BZI(I,J)= BZN(J)
		
	
400	   CONTINUE	
C-----------------------------------------------------------

C	VELOCITY FORM MOMENTUM
	
	
	
	DO J=1,INY

	GAMAI(I,J)=DSQRT(1.0D0+((PEXI(I,J))**2+(PEYI(I,J))**2+
     &    (PEZI(I,J))**2))

	VELXI(I,J)= PEXI(I,J)/GAMAI(I,J)
	
	VELYI(I,J)= PEYI(I,J)/GAMAI(I,J)

	VELZI(I,J)= PEZI(I,J)/GAMAI(I,J)	

	ENDDO



C
C-------------------------------------------------------------
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
 800	      CONTINUE         ! CONTINUATION OF LOOP IN X-DIRECTION 
c========================================================================
C

!	call PHI2D(DENI,PHII)

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 900	 CONTINUE              !CONTINUATION OF HALF-FULL TIME STEP
C-------------------------------------------------------------------
C	
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	FINAL RESULT AFTER COMPLITION OF TIME STEP
C-----------------------------------------------------------------------
	   DO 901 I = 1, INX
	      DO 902 J = 1, INY
		DEN(I,J) = DENI(I,J)
		DENB(I,J) = DENBI(I,J)

		VELX(I,J) = VELXI(I,J)
		VELY(I,J) = VELYI(I,J)
		VELZ(I,J) = VELZI(I,J)
		
		PEX(I,J) = PEXI(I,J)
		PEY(I,J) = PEYI(I,J)
		PEZ(I,J) = PEZI(I,J)


		EX(I,J)= EXI(I,J)
		EY(I,J)= EYI(I,J)
		EZ(I,J)= EZI(I,J)

		BX(I,J)= BXI(I,J)
		BY(I,J)= BYI(I,J)
		BZ(I,J)= BZI(I,J)
		
 902	      CONTINUE
 901	   CONTINUE
C----------------------------------------------------------------------

	

	RETURN
	END SUBROUTINE GASD1
C======================================================================

