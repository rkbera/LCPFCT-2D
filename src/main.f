C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C**************PROGRAME FOR INVICID FLUID (EULAR EQUATION)**************
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	PROGRAM MAIN
	INCLUDE 'param.i'
	
C	GRID VARIABLES
	
	INTEGER NXP, NYP,NZP
	INTEGER NX, NY,NZ
	INTEGER  M, N,NP
	
	REAL*8 XCELL(NPTX), YCELL(NPTY)
	REAL*8 XINT(NXINT), YINT(NYINT)
	
	REAL*8 DX, DY,DZ
	REAL*8 LX, LY,LZ

	REAL*8 KXA(NPTX),KYA(NPTY),RAND(NPT),PP(NPT)
C------
	REAL*8 KX,KY,PI,SUM1,SUM2,SUM3,RAND2(NPT)
	
C----
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
C	COMPONENTS OF WORTICITY
	
	REAL*8 DEN(NPTX,NPTY),DENB(NPTX,NPTY),DENBS(NPTX,NPTY)
C------
C	COMPONENTS OF VELOCITY(V).
	
	REAL*8 VX, VY,DIVB(NPTX,NPTY)
	REAL*8 VELX(NPTX,NPTY), VELY(NPTX,NPTY),VELZ(NPTX,NPTY)
	REAL*8 VBLX(NPTX,NPTY), VBLY(NPTX,NPTY), VBLZ(NPTX,NPTY)
	REAL*8 VBSLX(NPTX,NPTY), VBSLY(NPTX,NPTY), VBSLZ(NPTX,NPTY)
	REAL*8 PHI(NPTX,NPTY),ELECY(NPTX,NPTY),ELECZ(NPTX,NPTY)
	REAL*8 ELEC(NPTX,NPTY),ELECX(NPTX,NPTY),GAMA(NPTX,NPTY)
	REAL*8 BX(NPTX,NPTY),BY(NPTX,NPTY),BZ(NPTX,NPTY),B(NPTX,NPTY)
	REAL*8 EF(NPTX,NPTY),KE(NPTX,NPTY)
	REAL*8 ENB(NPTX,NPTY),GAMAB(NPTX,NPTY)

	REAL*8 EFF(NPT),KEF(NPT),ENBF(NPT)
C------

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	TIME-EVOLUTION-VARIABLES/DO LOOP INDICES
	
	INTEGER I, J, K,TSTEP,P
C------
C	TIMESTEP RELATED PARAMETERS
	
	INTEGER INISTP, MAXSTP, JSTEP
C------
C	PRINTING RELATED PARAMETERS
	
	INTEGER LOUT, IPRINT,LOUTP
C------
C	PARAMETERS TO ADJUST TIME-STEP DT ACCORDING TO COURANT NUMBER
	
	REAL*8 VEL(NPTX,NPTY), COURNT, VMAX, DTNEW, DT, DTMAX, TIME
C------
C	LOGICAL VARIABLE TO RESTART THE CODE
	
	LOGICAL  RESTRT
C------

C----------------PARTICLE PARAMETERS
	INTEGER HX1,HY1,HX2,HY2

	REAL*8 XP0(NPT),YP0(NPT),XPN(NPT),YPN(NPT)
	REAL*8 VPX0(NPT),VPY0(NPT),VPXN(NPT),VPYN(NPT),VPZ0(NPT),VPZN(NPT)
	REAL*8 PPX0(NPT),PPXN(NPT),PPY0(NPT),PPYN(NPT),PPZ0(NPT),PPZN(NPT)
	REAL*8 GAMMP(NPT),EP(NPT)
	
	REAL*8 DELX(NPT),DELXP,XNEAR1,XNEAR2,SIGNX
	REAL*8 DELY(NPT),DELYP,YNEAR1,YNEAR2,SIGNY

	REAL*8 EPX,EPY,EPZ,BPX,BPY,BPZ
	real*8 DTH

C---------------------------------------------
      character(len=100) :: DIR, DIR_2, FILENAME
      character(len=20) :: timestamp
      
      CHARACTER*100 LINE
      CHARACTER*20 PARAM_NAME
      CHARACTER*20 PARAM_VALUE
      INTEGER IOS, EQ_POS

      INTEGER BEAM, TPS

      REAL*8 density_beam, velocity_beam, DENB0, VB0
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	COMMON-GLOBLE-BLOCK	
	COMMON / DATAP / NX,NY,NXP,NYP
	COMMON / VARIAB /VELX,VELY,DEN,DENB,VBLX,VBLY
	COMMON / GRID / XINT,YINT,XCELL,YCELL,DX,DY,LX,LY,PI
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        print*, 'Simulation starts....'

C=============================================================
C	INITIALIZING VARIABLES FROM INPUT FILE

        ! default values
	LX = 10.0
	LY = 10.0
	LOUT = 10
	IPRINT= 10
	COURNT = 0.1
        INISTP = 0
	MAXSTP = 100
	DTMAX = 1.0
	DT = 0.1
        

        BEAM = 0

        TPS = 0
	NP = 100
	LOUTP = 10

          DO  I = 1,NX
           DO J = 1,NY
                
                DEN(I,J) = 1.0D0        

                ELECX(I,J)=0.0D0
                ELECY(I,J)=0.0D0
                ELECZ(I,J)=0.0D0

                VELX(I,J)= 0.0D0        
                VELY(I,J)= 0.0D0
                VELZ(I,J)= 0.0D0

                BX(I,J)= 0.0D0
                BY(I,J)= 0.0D0
                BZ(I,J)= 0.0d0
         
                DENB(I,J) = 0.0D0
                VBLX(I,J)= 0.0D0
                VBLY(I,J)= 0.0D0
                VBLZ(I,J)= 0.0D0

           ENDDO
         ENDDO



C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	LOADING INPUT FROM FILE
	OPEN(UNIT=10, FILE='input_file', STATUS='OLD', ACTION='READ')

10    CONTINUE
      READ(10, '(A)', IOSTAT=IOS) LINE
      IF (IOS /= 0) GOTO 20

      CALL PARSE_LINE(LINE, PARAM_NAME, PARAM_VALUE, EQ_POS)

      IF (EQ_POS .GT. 0) THEN
          IF (PARAM_NAME .EQ. 'LX') THEN
              READ(PARAM_VALUE, *) LX
          ELSE IF (PARAM_NAME .EQ. 'LY') THEN
              READ(PARAM_VALUE, *) LY
          ELSE IF (PARAM_NAME .EQ. 'NX') THEN
              READ(PARAM_VALUE, *) NX
          ELSE IF (PARAM_NAME .EQ. 'NY') THEN
              READ(PARAM_VALUE, *) NY
          ELSE IF (PARAM_NAME .EQ. 'TIME') THEN
              READ(PARAM_VALUE, *) TIME
          ELSE IF (PARAM_NAME .EQ. 'COURNT') THEN
              READ(PARAM_VALUE, *) COURNT
          ELSE IF (PARAM_NAME .EQ. 'MAXSTP') THEN
              READ(PARAM_VALUE, *) MAXSTP
          ELSE IF (PARAM_NAME .EQ. 'INISTP') THEN
              READ(PARAM_VALUE, *) INISTP   
          ELSE IF (PARAM_NAME .EQ. 'DT') THEN
              READ(PARAM_VALUE, *) DT
          ELSE IF (PARAM_NAME .EQ. 'DTMAX') THEN
              READ(PARAM_VALUE, *) DTMAX
          ELSE IF (PARAM_NAME .EQ. 'LOUT') THEN
              READ(PARAM_VALUE, *) LOUT
          ELSE IF (PARAM_NAME .EQ. 'IPRINT') THEN
              READ(PARAM_VALUE, *) IPRINT

          ELSE IF (PARAM_NAME .EQ. 'BEAM') THEN
              READ(PARAM_VALUE, *) BEAM          
          ELSE IF (PARAM_NAME .EQ. 'density_beam') THEN
              READ(PARAM_VALUE, *) density_beam
          ELSE IF (PARAM_NAME .EQ. 'velocity_beam') THEN
              READ(PARAM_VALUE, *) velocity_beam

          ELSE IF (PARAM_NAME .EQ. 'TPS') THEN
              READ(PARAM_VALUE, *) TPS
          ELSE IF (PARAM_NAME .EQ. 'NP') THEN
              READ(PARAM_VALUE, *) NP
          ELSE IF (PARAM_NAME .EQ. 'LOUTP') THEN
              READ(PARAM_VALUE, *) LOUTP
          END IF
      END IF
      GOTO 10

20    CONTINUE
      CLOSE(10)

!------ check input values
      !PRINT *, 'BEAM =', BEAM
      !PRINT *, 'TPS =',  TPS
     
	!call exit(123)

C--------------------------------------------------------
	NX = NPTX
	NY = NPTY
        
	NXP = NX+1
	NYP = NY+1	

	PI = 4.0D0*DATAN(1.0D0)


!C------------INITIALIZE TEST PARTICLE PARAMETER for TPS ----------

       IF (TPS .EQ. 1) THEN

       ! Generate random numbers for distributing test particle randomly
       call GENERATE_RANDOM(NP, RAND)
       call GENERATE_RANDOM_2(NP, RAND2)

C-------------------------------------
!! Initial position of the particle

	DO P = 1, NP
!! Initial position of the particle

        XP0(P)= LX*RAND(P);
        YP0(P)= LY*RAND2(P); 

	VPX0(P)=0.0d0; 
	VPY0(P)=0.0d0; 
	VPZ0(P)=0.0d0

	GAMMP(P)=1.0d0/dsqrt(1.0d0-VPX0(P)**2-VPY0(P)**2-VPZ0(P)**2)

!------------MOMENTUM
	PPX0(P)=GAMMP(P)*VPX0(P)
	PPY0(P)=GAMMP(P)*VPY0(P)
	PPZ0(P)=GAMMP(P)*VPZ0(P)
	
	EP(P)=GAMMP(P)  !ENERGY

	ENDDO
        
       END IF

!-------------------
	EPX=0.0d0;
!------------------	
	
	DX = (LX)/DFLOAT(NX)
	
	DY = (LY)/DFLOAT(NY)
	
	
	
	RESTRT = .FALSE.
	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	SETTING CELL-INTERFACE LOCATIONS
	
	DO 1 I = 1, NXP
		XINT(I) = DFLOAT(I-1)*DX
 1	CONTINUE
	
	DO 2 J = 1, NYP
		YINT(J) = -LY/2.0D0+ DFLOAT(J-1)*DY
 2	CONTINUE

	
C-------
C	SETTING CELL-CENTRE LOCATIONS
	
	DO 4 I = 1, NX
		XCELL(I) = XINT(I) + DX/2.0D0
 4	ENDDO
	
	DO 5 J = 1, NY
		YCELL(J) = YINT(J) + DY/2.0D0
 5	ENDDO

	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	CALLING SUBROUTINE INDATA FOR INITIAL PROFILES 
	
	CALL PLASMA_INIPROF(DEN,VELX,VELY,VELZ,ELECX,ELECY,ELECZ,BX,BY,BZ)
      
        IF (BEAM .EQ. 1) THEN
        DENB0 = density_beam
        VB0  = velocity_beam
        CALL BEAM_INIPROF(DENB,VBLX,VBLY,VBLZ,DENB0,VB0)
        END IF
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	DIVB(1:NX,1)=0.0D0; DIVB(1,1:NY)=0.0D0;DIVB(NX,1:NY)=0.0D0;
	DIVB(1:NX,NY)=0.0D0;

	DO I=2,NX-1
	DO J=2,NY-1

	DIVB(I,J)=((BX(I+1,J)-BX(I-1,J))/2*DX)+(BY(I,J+1)-BY(I,J-1))/2*DY

	ENDDO
	ENDDO


	DO I=1,NX
	DO J=1,NY


	GAMAB(I,J)=1.0D0/DSQRT(1.0D0 -(VBLX(I,J))**2-
     &	(VBLY(I,J))**2-(VBLZ(I,J))**2)
	
	ENB(I,J)=DENB(I,J)*(GAMAB(I,J)-1.0D0)

	GAMA(I,J)=1.0D0/DSQRT(1.0D0 -(VELX(I,J))**2-
     &	(VELY(I,J))**2-(VELZ(I,J))**2)
	
	KE(I,J)=DEN(I,J)*(GAMA(I,J)-1.0D0)
	
	EF(I,J)=((ELECX(I,J))**2) +((ELECY(I,J))**2) +((ELECZ(I,J))**2)
     &   +((BX(I,J))**2) +((BY(I,J))**2) +((BZ(I,J))**2) 

	ENDDO
	ENDDO

c--------------------------------------
	!OPEN(UNIT=55,file='energy.dat',STATUS='UNKNOWN')

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         print*, 'Initialization successful!'
         print*, 'Entering Time loop ...'

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	BEGINING OF TIME STEPS
	
	TIME =0.0D0

	DO 13 TSTEP = INISTP, MAXSTP
!--------------------------------------------
         print*, 'Timestep =', TSTEP-1, ', dt =', DT, ', Time =', TIME
C---------------------------------------------------------------------
C		ADJUSTING DT TO SATISFY COURANT CONDITION
	
		VMAX = 0.0D0
	   DO 15 J = 1, NY
	      DO 16 I = 1, NX
		VEL(I,J)=DSQRT((VELX(I,J))**2 + (VELY(I,J))**2)
		VMAX = MAX(VELX(I,J),VMAX)
 16	      CONTINUE
 15	   CONTINUE
	
	
	DTNEW = COURNT*MIN(DX,DY)/(VMAX)
	
!	DT= MIN(DTNEW,DT,DTMAX) 

	DT= (COURNT*DX)/DSQRT(2.0D0)
C--------------------------------------------------------------------
C	PRINTING THE RESULT AS REQUIRED
	
	JSTEP = TSTEP-1
	
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C-------- ENERGY CALCULATION

	SUM1=0.0D0
	SUM2=0.0D0
	SUM3=0.0D0

	DO I=1,NX
	DO J=1,NY

	SUM1=SUM1+KE(I,J)

	SUM2=SUM2+EF(I,J)

	SUM3=SUM3+ENB(I,J)
	ENDDO
	ENDDO

	KEF(TSTEP)=SUM1
	EFF(TSTEP)=SUM2
	ENBF(TSTEP)=SUM3
		

!	WRITE(55,*) TIME,EFF(TSTEP),KEF(TSTEP),ENBF(TSTEP)


C===========================================
C       creating output directory to save data
        ! Define the directory
	DIR = 'PLASMA_DATA'
        
        ! Create the directory using a system call
	call system('mkdir -p ' // trim(DIR))

        IF (TPS .EQ. 1) THEN
	 ! Define the directory
        DIR_2 = 'TEST_PARTICLE_DATA'

        ! Create the directory using a system call
        call system('mkdir -p ' // trim(DIR_2))
        
        END IF
c------------------------------------
	
	IF ( (JSTEP/IPRINT)*IPRINT .EQ. JSTEP ) THEN
C-------
        print*, 'writing plasma data...'

C	PRINTING THE VALUES OF VELX, VELY, ND AT DIFFERENT TIMES
	! Generate a timestamp or step number for the filename
	write(timestamp, '(I4.4)') JSTEP ! Format step number with leading zeros

        ! Construct the filename with the directory and timestamp
        FILENAME = trim(DIR)// '/plasma_data_' //trim(timestamp)//'.dat'

        ! Open the file for writing
	open(unit=10,file=trim(FILENAME),status='replace',action='write')

        
	DO 18 J = 1, NY
	DO 19 I = 1, NX

	WRITE(10,*) XCELL(I), YCELL(J), DEN(I,J),DENB(I,J),ELECX(I,J)

 19     CONTINUE
 18     CONTINUE

        ! Close the file
	close(10)


        IF (TPS .EQ. 1) THEN

        print*, 'writing particle data...'
        
	! Construct the filename with the directory and timestamp
        FILENAME = trim(DIR_2)// '/tp_data_' //trim(timestamp)//'.dat'

        ! Open the file for writing
	open(unit=20,file=trim(FILENAME),status='replace',action='write')
	!CALL PRINT4 ('tp_data',JSTEP,LOUTP)

	DO P= 1, NP     
	WRITE(20,*)  XP0(P),YP0(P), VPX0(P), VPY0(P), EP(P)
	ENDDO

	 ! Close the file
        close(20)

	END IF

	
	ENDIF
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	CALLING GASDYNAMICS/AFTER TIMESTEP DT.
	
	TIME = TIME + DT
C**************************************************************

	CALL  GASD1(DT,VELX,VELY,VELZ,DEN,VBLX,VBLY,VBLZ,DENB,
     &   ELECX,ELECY,ELECZ,BX,BY,BZ,TSTEP) 

C---------------------------------------------------
	DIVB(1:NX,1)=0.0D0; DIVB(1,1:NY)=0.0D0;DIVB(NX,1:NY)=0.0D0;
	DIVB(1:NX,NY)=0.0D0;

	DO I=2,NX-1
	DO J=2,NY-1

	DIVB(I,J)=((BX(I+1,J)-BX(I-1,J))/2*DX)+(BY(I,J+1)-BY(I,J-1))/2*DY
	
	ENDDO
	ENDDO

	DO I=1,NX
	DO J=1,NY


	GAMAB(I,J)=1.0D0/DSQRT(1.0D0 -(VBLX(I,J))**2-
     &	(VBLY(I,J))**2-(VBLZ(I,J))**2)
	
	ENB(I,J)=DENB(I,J)*(GAMAB(I,J)-1.0D0)

	GAMA(I,J)=1.0D0/DSQRT(1.0D0 -(VELX(I,J))**2-
     &	(VELY(I,J))**2-(VELZ(I,J))**2)
	
	KE(I,J)=DEN(I,J)*(GAMA(I,J)-1.0D0)
	
	EF(I,J)=((ELECX(I,J))**2) +((ELECY(I,J))**2) +((ELECZ(I,J))**2) 
     &   +((BX(I,J))**2) +((BY(I,J))**2) +((BZ(I,J))**2) 

	ENDDO
	ENDDO

C***********************************************************************
C---- TEST PARTICLE SIMULATION -- ADVANCE PATICLES POS and VEL-------
       
        IF (TPS .EQ. 1) THEN

	DO 1000 P=1,NP

!---------------------------------------------
!   finding nearest cell for the particles 

!-------along x
	DO I=1,NX	

	DELX(i)=ABS(XP0(P)-XCELL(I))

	ENDDO

	DELXP=MINVAL(DELX(1:NX))


	XNEAR1=XP0(P) - DELXP ! nearest grid point

	HX1=INT(ABS(XNEAR1)/DX)

!-------- along y
	DO J=1,NY

	DELY(J)=ABS(YP0(P)-YCELL(J))

	ENDDO

	DELYP=minval(DELY(1:NY))


	YNEAR1=YP0(P) - DELYP ! nearest grid point

	HY1=INT(abs(YNEAR1)/DY)

c--------finding particle postion right or left to the minimum grid
	
	IF (DELXP .ge. 0.0d0) THEN

	SIGNX=1  ! right to the grid

	ELSE

	SIGNX=0  !left to the grid

	ENDIF


	IF (SIGNX .EQ. 1) THEN

	XNEAR2= XNEAR1 + DX

	HX2=HX1+1
	
	ELSE

	XNEAR2=XNEAR1-DX

	HX2=HX1-1

	ENDIF

c--------finding particle postion right or left to the minimum grid
	
	IF (DELYP .ge. 0.0d0) THEN

	SIGNY=1  ! right to the grid

	ELSE

	SIGNY=0  !left to the grid

	ENDIF

	IF (SIGNY .EQ. 1) THEN

	YNEAR2=YNEAR1+ DY

	HY2=HY1+1
	

	ELSE

	YNEAR2=YNEAR1-DY

	HY2=HY1-1

	ENDIF

c-----------------------
! field interpolation at the particle position
!---------------------------
	EPX=ELECX(HX1,HY1)+ELECX(HX1,HY2)+ELECX(HX2,HY1)+ELECX(HX2,HY2)
	EPY=ELECY(HX1,HY1)+ELECY(HX1,HY2)+ELECY(HX2,HY1)+ELECY(HX2,HY2)
	EPZ=ELECZ(HX1,HY1)+ELECZ(HX1,HY2)+ELECZ(HX2,HY1)+ELECZ(HX2,HY2)

	EPX=EPX/4.0D0;
	EPY=EPY/4.0D0;
	EPZ=EPZ/4.0D0;
	
	BPX=(BX(HX1,HY1)+BX(HX1,HY2)+BX(HX2,HY1)+BX(HX2,HY2))/4.0D0
	BPY=(BY(HX1,HY1)+BY(HX1,HY2)+BY(HX2,HY1)+BY(HX2,HY2))/4.0D0
	BPZ=(BZ(HX1,HY1)+BZ(HX1,HY2)+BZ(HX2,HY1)+BZ(HX2,HY2))/4.0D0

C------------- UPGRADE PARTICLE POSITION
	
  
	DTH = 0.5d0*DT

	PPXN(P)=PPX0(P) -DTH*(EPX+(VPY0(P)*BPZ) -(VPZ0(P)*BPY))

	PPYN(P)=PPY0(P) -DTH*(EPY+(VPZ0(P)*BPX) -(VPX0(P)*BPZ))

	GAMMP(P)=DSQRT(1.0d0+(PPXN(P)**2)+(PPYN(P)**2))

	VPXN(P)=PPXN(P)/GAMMP(P)
	VPYN(P)=PPYN(P)/GAMMP(P)


	XPN(P)=XP0(P)+DT*VPXN(P)
	YPN(P)=YP0(P)+DT*VPYN(P)

!----upgrade particle 

	XP0(P)=XPN(P)
	VPX0(P)=VPXN(P)
	PPX0(P)=PPXN(P)

	YP0(P)=YPN(P)
	VPY0(P)=VPYN(P)
	PPY0(P)=PPYN(P)

	EP(P)=GAMMP(P)

1000	CONTINUE
	
        END IF
!------------------------------------------

	
 13	CONTINUE                  !CONTINUATION OF TIME STEP

	print*, "Simulation ended successfully!"
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	END PROGRAM
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
