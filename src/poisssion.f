c-----------------------------------------------------------------
c	for 1D potential	
	SUBROUTINE phi_e(k1,kn,dx,rhon,e)



	implicit real*8( a-h, o-z )

	parameter( npt=3000 )

	real*8 phi(npt),e(npt),r(npt),phi0,atd(npt),btd(npt),

     &	   ctd(npt),phiav,ph(npt),rhon(kn),pot(npt)



	do 2 i = 1, kn - 1

	   r(i) = 0.0d0

2	continue

	phi0 = 0.0d0

	phi(kn) = phi0

	do 3 i1 = 1, kn - 1

	   btd(i1) = -2.0d0

3	continue     

	do 4 i2 = 1, kn - 2

	   atd(i2+1) = 1.0d0

	   ctd(i2)   = 1.0d0

4	continue

	do 5 i3 = 1, kn - 1

	      r(i3) =  r(i3) +  (rhon(i3) - 1.0d0)*dx**2

5	continue

	r(1)    = r(1) - phi0
	r(kn-1) = r(kn-1) - phi0
	r(kn) = 0.0
!-------------------------------------------------------------

	call tridag( atd, btd, ctd, r, ph, kn-1 )

!-------------------------------------------------------------

	do 8 i4 = k1, kn-1

	   phi(i4) = ph(i4)

8	continue



!-------------Subtract out the average-------------------------
 

	phiav = 0.0d0

	do 9 i5 = k1, kn 

	   phiav = phi(i5) + phiav

9	continue

	phiav = phiav / dfloat(kn-k1+1)

	do 10 i6 = k1, kn

	   phi(i6) = phi(i6) - phiav
           pot(i6)=phi(i6)
          
        
10	continue


!-------------Calculate electric field---------------------------



	do 11 i7 = k1+1, kn-1

	   e(i7) = -( phi(i7+1) - phi(i7-1) ) / (2.0d0*dx)

11	continue

	e(k1) = -( phi(k1+1) - phi(kn) ) / (2.0d0*dx)

	e(kn) = -( phi(k1) - phi(kn-1) ) / (2.0d0*dx)



	return

	end



!*******************************************************************



      SUBROUTINE TRIDAG(A,B,C,R,U,N)

	implicit real*8 ( a-h, o-z )

      PARAMETER (NMAX=3000)

      DIMENSION GAM(NMAX),A(N),B(N),C(N),R(N),U(N)

!      IF(B(1).EQ.0.0D0)PAUSE

      BET=B(1)

      U(1)=R(1)/BET

      DO 11 J=2,N

        GAM(J)=C(J-1)/BET

        BET=B(J)-A(J)*GAM(J)

!       IF(BET.EQ.0.0D0)PAUSE

        U(J)=(R(J)-A(J)*U(J-1))/BET

11    CONTINUE

      DO 12 J=N-1,1,-1

        U(J)=U(J)-GAM(J+1)*U(J+1)

12    CONTINUE

      RETURN

      END       

c=================================================================
C************** POISSION EQUTION 2D *************************** |


	subroutine PHI2D(DEN,ELECX,ELECY)

C	GRID VARIABLES 
	INCLUDE 'param.i'
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
	
	REAL*8 DEN(NPTX,NPTY),DENB(NPTX,NPTY)
	REAL*8 PHI(NPTX,NPTY),RHO(NPTX,NPTY)
	REAL*8 ELECY(NPTX,NPTY),PI
	REAL*8 ELEC(NPTX,NPTY),ELECX(NPTX,NPTY),ELECZ(NPTX,NPTY)
	REAL*8 BX(NPTX,NPTY), BY(NPTX,NPTY), BZ(NPTX,NPTY)
	REAL*8 DENX(NPTX),EXN(NPTX),EYN(NPTY),DENY(NPTY)
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
	
	XLEFT  = XCELL(1)
	XRIGHT = XCELL(NPTX)+DX
	YBOTOM = YCELL(1)
	YTOP   = YCELL(NPTY)+DY
C~~~~~~~~~~~~~~~~~~~~~~~FOR_DIPOLE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	DO 110 I = 1,NX
	   DO 120 J = 1,NY
	    RHO(I,J)=DEN(I,J)-1.0D0
 120 	CONTINUE
 110	CONTINUE
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	TO CALCULATE VELOCITY USING HELMHOLTZ SOLVER

	CALL HELM(RHO,PHI,XLEFT,XRIGHT,YBOTOM,YTOP)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	

	DO I=2,NX-1
	DO J=2,NY-1
	
	ELECX(I,J)= -(PHI(I+1,J)-PHI(I-1,J))/(2*DX)
	
	ENDDO
	ENDDO

	DO J=2,NY-1
	DO I=2,NX-1
	
	ELECY(I,J)= -(PHI(I,J+1)-PHI(I,J-1))/(2*DY)
	
	ENDDO
	ENDDO
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	DO J=1,NY-1
	ELECX(NX,J)=ELECX(NX-1,J)+(ELECX(NX-1,J)-ELECX(NX-2,J))
	ELECY(NX,J)=ELECY(NX-1,J)+(ELECY(NX-1,J)-ELECY(NX-2,J))

	ELECX(1,J)=(ELECX(NX,J)+ELECX(2,J))/2.0D0
	ELECY(1,J)=(ELECY(NX,J)+ELECY(2,J))/2.0D0
	ENDDO
	

	DO I=1,NX-1
	ELECX(I,NY)=ELECX(I,NY-1)+(ELECX(I,NY-1)-ELECX(I,NY-2))
	ELECY(I,NY)=ELECY(I,NY-1)+(ELECY(I,NY-1)-ELECY(I,NY-2))

	ELECX(I,1)=(ELECX(I,NY)+ELECX(I,2))/2.0D0
	ELECY(I,1)=(ELECY(I,NY)+ELECY(I,2))/2.0D0
	ENDDO
	
	ELECX(NX,NY)= ELECX(NX-1,NY-1)
	ELECY(NX,NY)= ELECY(NX-1,NY-1)
C====================================================================
	
	
! 	DO J=1,NY-1
!	

!	ELECX(NX,J)=(ELECX(NX-1,J)+ELECX(1,J))/(2)
!	ELECY(NX,J)=(ELECY(NX-1,J)+ELECY(1,J))/(2)
!	
!	ENDDO

!	DO I=1,NX-1
!	

!	ELECX(I,NY)=(ELECX(I,NY-1)+ELECX(I,1))/(2)
!	ELECY(I,NY)=(ELECY(I,NY-1)+ELECY(I,1))/(2)
!	
!	ENDDO

	



	end subroutine PHI2D
  

