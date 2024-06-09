c  Subroutines for LCPFCT package.
c  Some or all of these routines must be called by the calling program.
c  Last modified: 7/2/94.

	SUBROUTINE gasdyn (k1,kn,bc1,bcn,dt)

c++++++++++++++++++++++++++++++++++++++
c  Description: This routine integrates the gasdynamic equations using
c  the momentum component RVRN as the direction of integration and the
c  momentum RVTN as the transverse direction. In 2D models, the
c  two directions of integration are chosen by exchanging RVRN and
c  RVTN in Common.

c  Arguments:
c  (k1, kn) -- I/P integer: Index of the integration's (first, last) cell
c  (bc1, bcn) : I/P integer: indicates BC on integration at (k1, kn).
c  (dt) : I/P real: timestep for integration of this step.
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,k1,k1p,bc1,bcn,k,kn,knp,it
	parameter (npt=200000)
	logical pbc
	real*8 sbc1,srv1,sbcn,srvn,vrho1,vrhon,vrvr1,vrvrn,vrvt1,
     &	vrvtn,verg1,vergn,mpint(npt),vel(npt),unit(npt),zero(npt),
     &	rhoo(npt),rvro(npt),rvto(npt),ergo(npt),vint(npt),pre(npt),
     &	mpvint(npt),dtsub,dt,relax

	real*8 rho_in,pre_in,vel_in,gamma0,rhoamb,preamb,velamb,gammam,
     &	rhon(npt),rvrn(npt),rvtn(npt),ergn(npt)

	common /arrays/ rhon,rvrn,rvtn,ergn,relax,rho_in,pre_in,vel_in,
     &	gamma0,rhoamb,preamb,velamb,gammam
     
	data unit/npt*1.0d0/, zero/npt*0.0d0/     

c  Prepare for time integration. Index K is either I or J depending
c  on the definitions of RVRN and RVTN. Copies of the physical variable
c  are needed to recover values for the whole step integration.

	knp = kn + 1
	k1p = k1 + 1
	pbc = .false.
	if (bc1.eq.3 .or. bcn.eq.3) pbc = .true.
	do 50 k = k1,kn
		rhoo(k) = rhon(k)
		rvro(k) = rvrn(k)
		rvto(k) = rvtn(k)
50		ergo(k) = ergn(k)

c  Integrate first the half step, then the whole step:
	do 500 it = 1,2
	dtsub = 0.5d0 * dt * dfloat(it)
	do 100 k = k1,kn
		vel(k) = rvrn(k) / rhon(k)
100		pre(k) = gammam * (ergn(k) - 0.5d0 *
     &			 (rvrn(k)**2 + rvtn(k)**2)/rhon(k))

c  Calc. the interface velocities and pressures as weighted values
c  of the cell-centered values computed just above:
	do 200 k = k1+1,kn
		mpvint(k) = 1.0d0 / (rhon(k) + rhon(k-1))
		vint(k) = (vel(k)*rhon(k-1) + vel(k-1)*rhon(k)) * 
     &			  mpvint(k)
		mpint(k) = -(pre(k)*rhon(k-1) + pre(k-1)*rhon(k)) *
     &			   mpvint(k)
200		mpvint(k) = -(pre(k)*vel(k)*rhon(k-1) + pre(k-1) *
     &			    vel(k-1)*rhon(k)) * mpvint(k)

c  The unweighted interface averages can be computed as follows:
c	vint(k) = 0.5d0 * (vel(k) + vel(k-1))
c	mpint(k) = -0.5d0 * (pre(k) + pre(k-1))
c200	mpvint(k) = mpint(k) * vint(k)

c  Call the FCT utility routines and set the BC. Other BC could be
c  added for inflow, outflow, etc.
c  (BC1, BCN) = 1  => ideal solid wall or axis boundary condition
c  (BC1, BCN) = 2  => an extrapolative outflow BC
c  (BC1, BCN) = 3  => periodic BC
c  (BC1, BCN) = 4  => specified boundary values (e.g. shock tube problem)

c+++++  Assign BC at "k1"  ++++++++++++
	go to (310,320,330,340), bc1

310	vint(k1) = 0.0d0
	mpint(k1) = -pre(k1)
	mpvint(k1) = 0.0d0
	go to 350

320	vint(k1) = vel(k1) * (1.0d0 - relax)
	mpint(k1) = -pre(k1) * (1.0d0-relax) - relax*pre_in
	mpvint(k1) = mpint(k1) * vint(k1)
	go to 350

330	mpvint(k1) = 1.0d0 / (rhon(k1) + rhon(kn))
	vint(k1) = mpvint(k1) * (vel(k1)*rhon(kn) + vel(kn)*rhon(k1))
	mpint(k1) = -mpvint(k1) * (pre(k1)*rhon(kn) + pre(kn)*rhon(k1))
	mpvint(k1) = -mpvint(k1) * (pre(k1)*vel(k1)*rhon(kn) + 
     &		      pre(kn)*vel(kn)*rhon(k1))
	go to 350

340	vint(k1) = vel_in
	mpint(k1) = -pre_in
	mpvint(k1) = -pre_in * vel_in

c+++++  Assign BC at "kn"  ++++++++++++
350	go to (410,420,430,440), bcn

410	vint(knp) = 0.0d0
	mpint(knp) = -pre(kn)
	mpvint(knp) = 0.0d0
	go to 450

420	vint(knp) = vel(kn) * (1.0d0 - relax)
	mpint(knp) = -pre(kn) * (1.0d0-relax) - relax * preamb
	mpvint(knp) = mpint(knp) * vint(knp)
	go to 450

430	vint(knp) = vint(k1)
	mpint(knp) = mpint(k1)
	mpvint(knp) = mpvint(k1)
	go to 450

440	vint(knp) = velamb
	mpint(knp) = -preamb
	mpvint(knp) = -preamb * velamb

450	continue
c++++++++++++++++++++++++++++++++++++++

c  The velocity-dependent FCT coefficients are set and the BC calc. 
c  are completed. Here the periodic BC require no action as (S)lope
c  and (V)alue boundary value specifiers are ignored in LCPFCT when 
c  PBC = .true.

	CALL velocity (vint,k1,knp,dtsub)

	go to (510,520,550,540), bc1

c++++++++
510	CALL zeroflux (k1)
	sbc1 = 1.0d0
	srv1 = -1.0d0
	vrho1 = 0.0d0
	vrvr1 = 0.0d0
	vrvt1 = 0.0d0
	verg1 = 0.0d0
	go to 550

520	CALL zerodiff (k1)
	sbc1 = 1.0d0 - relax
	srv1 = 1.0d0 - relax
	vrho1 = relax * rho_in
	vrvr1 = 0.0d0
	vrvt1 = 0.0d0
	verg1 = relax * pre_in / gammam
	go to 550
	
540	sbc1 = 0.0d0
	srv1 = 0.0d0
	vrho1 = rho_in
	vrvr1 = rho_in * vel_in
	vrvt1 = 0.0d0
	verg1 = pre_in / gammam + 0.5d0 * rho_in * vel_in**2

c++++++++
550	go to (610,620,650,640), bcn

610	CALL zeroflux (knp)
	sbcn = 1.0d0
	srvn = -1.0d0
	vrhon = 0.0d0
	vrvrn = 0.0d0
	vrvtn = 0.0d0
	vergn = 0.0d0
	go to 650

620	CALL zerodiff (knp)
	sbcn = 1.0d0 - relax
	srvn = 1.0d0 - relax
	vrhon = relax * rhoamb
	vrvrn = 0.0d0
	vrvtn = 0.0d0
	vergn = relax * preamb / gammam
	go to 650

640	sbcn = 0.0d0
	srvn = 0.0d0
	vrhon = rhoamb
	vrvrn = rhoamb * velamb
	vrvtn = 0.0d0
	vergn = preamb / gammam + 0.5d0 * rhoamb * velamb**2

650	continue

c+++++  Integrate the continity equations using LCPFCT  +++++
	CALL lcpfct (rhoo,rhon,k1,kn,sbc1,vrho1,sbcn,vrhon,pbc)

	CALL sources (k1,kn,dtsub,5,unit,mpint,mpint(k1),mpint(knp))

	CALL lcpfct (rvro,rvrn,k1,kn,srv1,vrvr1,srvn,vrvrn,pbc)

	CALL lcpfct (rvto,rvtn,k1,kn,sbc1,vrvt1,sbcn,vrvtn,pbc)

	CALL sources (k1,kn,dtsub,4,unit,mpvint,mpvint(k1),mpvint(knp))

	CALL lcpfct (ergo,ergn,k1,kn,sbc1,verg1,sbcn,vergn,pbc)

500	continue   !  End of halfstep-wholestep loop

	return
	END

c********************************************************************

	SUBROUTINE LCPFCT (rhoo,rhon,i1,in,srho1,vrho1,srhon,vrhon,pbc)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Originated: J.P. Boris        Code 4400, NRL,   Feb. 1987
c  Modified: Lab. for Computational Physics and Fluid Dynamics
c  Contact:  J.P. Boris, J.H. Gardner, A.M. Landsberg, E.S. Oran

c  DESCRIPTION: This routine solves generalized continuity equations
c  of the form: d(RHO)/dt = -div(rho*v) + SOURCES   in the users's 
c  choice of Cartesian, cylindrical or spherical coordinate systems.
c  A facility is included to allow definition of other coordinates. 
c  The grid can be Eulerian, sliding rezone, or Lagrangian, and can 
c  be arbitrarily spaced. The algorithm is a low phase-error FCT 
c  algorithm, vectorized and optimized for a combination of speed and 
c  flexibility. 
c  A complete description appears in the NRL Memorandum Report: 
c  "LCPFCT - A Flux-Corrected Transport Algorithm for Solving 
c  Generalized Continuity Equations", NRL/MR/6410-93-7192 (April 1993).


c  Arguments:
c  (rhoo) -- I/P real array: grid point densities at start of step
c  (rhon) -- O/P real array: grid point densities at end of step 
c  (i1, in) -- I/P integer: (first, last) grid point of integration
c  (srho1) -- I/P real array: boundary guard cell factor at cell (i1-1)
c  (vrho1) -- I/P real array: boundary value added to guard cell (i1-1)
c  (srhon) -- I/P real array: boundary guard cell factor at cell (in+1)
c  (vrhon) -- I/P real array: boundary value added to guard cell (in+1)
c  (pbc) -- logical: periodic boundaries if PBC = .true.

c  In this routine, the last interface at RADHN(INP) is the outer 
c  boundary of the last cell indexed "in". The first interface at 
c  RADHN(i1) is the outer boundary of the integration domain before 
c  the first cell indexed "i1".

c  Language and limitations:
c  LCPFCT is a package of FORTRAN 77 routines written in single 
c  precision (64 bits CRAY). The parameter "npt" is used to establish 
c  the internal FCT array dimensions at the maximum size expected.  
c  Thus npt=200000 means that continuity equations for systems up to 200 
c  cells long in one direction can be integrated. Underflows can occur
c  when the function being transported has a region of zeros. The 
c  calculations misconserve by one or two bits per cycle. Relative phase
c  and amplitude errors (for smooth functions) are typically a few 
c  percent for characteristic lengths of 1-2 cells (wavelengths of 
c  order 10 cells). The jump conditions for shocks are generally 
c  accurate to better than 1 percent. COMMON blocks are used to transmit
c  all data between routines in the package.

c  AUXILIARY ROUTINES:
c  CNVFCT, CONSERVE, COPYGRID, MAKEGRID, NEW_GRID, RESIDIFF, SET_GRID,
c  SOURCES, VELOCITY, ZERODIFF, ZEROFLUX.
c  The detailed documentation report provided (or the listing below)
c  explains the definitions and use of the arguments to these other
c  routines. These routines are not called from LCPFPT itself but are
c  controlled by calls from the user.  Routines MAKEGRID, VELOCITY and 
c  SOURCES must first be called to set the grid geometry, velocity-
c  dependent flux and diffusion coefficients, and external source 
c  arrays used by LCPFCT. The other routines may be called to perform
c  other functions such as to modify boundary conditions, to perform
c  special grid operations, or compute conservation sums.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i1,in,i1p,inp,i
	real*8 bignum,srho1,vrho1,srhon,vrhon,rho1m,rhonp,rhot1m,rhotnp,
     &	rhotd1m,rhotdnp
	logical pbc
	parameter(npt=200000, bignum=1.0d38)
c  (bignum) = machine-dependent largest number -- set by the user !!

	real*8 rhoo(npt),rhon(npt)

c  /FCT_SCRH/ holds scratch arrays for use by LCPFCT & CNVFCT:
	real*8 scrh(npt),scr1(npt),diff(npt),flxh(npt),fabs(npt),
     &	fsgn(npt),term(npt),terp(npt),lnrhot(npt),lorhot(npt),
     &	rhot(npt),rhotd(npt)
	common /fct_scrh/ scrh,scr1,diff,flxh,fabs,fsgn,term,terp,
     &	lnrhot,lorhot,rhot,rhotd

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo


c  /FCT_VELO/ holds velocity-dependent flux coefficients:
	real*8 hadudth(npt),nulh(npt),mulh(npt),epsh(npt),vdtodr(npt)
	common /fct_velo/ hadudth,nulh,mulh,epsh,vdtodr

c  /FCT_MISC/ holds the source array and diffusion coefficient:
	real*8 source(npt),diff1
	common /fct_misc/ source,diff1

c++++++++++++++++++++++++++++++++++++++
	i1p = i1 + 1
	inp = in + 1

c+++++  Calc. convective and diffusive fluxes  +++++
	if (pbc) then
		rho1m = rhoo(in)
		rhonp = rhoo(i1)
	else
		rho1m = srho1 * rhoo(i1) + vrho1
		rhonp = srhon * rhoo(in) + vrhon
	endif

	diff(i1) = nulh(i1) * (rhoo(i1) - rho1m)
	flxh(i1) = hadudth(i1) * (rhoo(i1) + rho1m)
	do 1 i = i1p,in
	flxh(i) = hadudth(i) * (rhoo(i) + rhoo(i-1))
1	diff(i) = nulh(i) * (rhoo(i) - rhoo(i-1))
	diff(inp) = nulh(inp) * (rhonp - rhoo(in))
	flxh(inp) = hadudth(inp) * (rhonp + rhoo(in))

c  Calc. "lorhot", the transported mass elements, and "lnrhot", the
c  transported & diffused mass elements:
	do 2 i = i1,in
	lorhot(i) = lo(i)*rhoo(i) + source(i) + (flxh(i)-flxh(i+1))
	lnrhot(i) = lorhot(i) + (diff(i+1) - diff(i))
	rhot(i) = lorhot(i) * rlo(i)
2	rhotd(i) = lnrhot(i) * rln(i)

c  Evaluate the boundary conditions for (rhot, rhotd):
	if (pbc) then
		rhot1m = rhot(in)
		rhotnp = rhot(i1)
		rhotd1m = rhotd(in)
		rhotdnp = rhotd(i1)
	else
		rhot1m = srho1 * rhot(i1) + vrho1
		rhotnp = srhon * rhot(in) + vrhon
		rhotd1m = srho1 * rhotd(i1) + vrho1
		rhotdnp = srhon * rhotd(in) + vrhon
	endif

c++++++++++++++++++++++++++++++++++++++
c  Calc. the transported antidiffusive fluxes and transported and
c  diffused density differences:
	flxh(i1) = mulh(i1) * (rhot(i1) - rhot1m)
	diff(i1) = rhotd(i1) - rhotd1m
	fabs(i1) = dabs(flxh(i1))
	fsgn(i1) = sign(diff1,diff(i1))

	do 3 i = i1p,in
	flxh(i) = mulh(i) * (rhot(i) - rhot(i-1))
3	diff(i) = rhotd(i) - rhotd(i-1)

	flxh(inp) = mulh(inp) * (rhotnp - rhot(in))
	diff(inp) = rhotdnp - rhotd(in)
c++++++++++++++++++++++++++++++++++++++

c  Calc. the magnitude & sign of the antidiffusive flux followed
c  by the flux-limiting changes on the right and left:
	do 4 i = i1,in
	fabs(i+1) = dabs(flxh(i+1))
	fsgn(i+1) = sign(diff1,diff(i+1))
	term(i+1) = fsgn(i+1) * ln(i) * diff(i)
4	terp(i) = fsgn(i) * ln(i) * diff(i+1)

	if (pbc) then
		terp(inp) = terp(i1)
		term(i1) = term(inp)
	else
		terp(inp) = bignum
		term(i1) = bignum
	endif

c  Correct the transported fluxes completely and then calculate the new
c  Flux-Corrected Transport densities:
	flxh(i1) = fsgn(i1) * dmax1(0.0d0,dmin1(term(i1),fabs(i1),
     &						terp(i1)))
	do 5 i = i1,in
	flxh(i+1) = fsgn(i+1) * dmax1(0.0d0,dmin1(term(i+1),fabs(i+1),
     &						  terp(i+1)))
	rhon(i) = rln(i) * (lnrhot(i) + (flxh(i) - flxh(i+1)))
5	source(i) = 0.0d0

	return
	END

c********************************************************************

	SUBROUTINE makegrid (radho,radhn,i1,inp,alpha)

c++++++++++++++++++++++++++++++++++++++
c  Description: Initializes geometry variables and coefficients. It 
c  should be called first to initialize the grid. The grid must be 
c  defined for all of the grid interfaces from "I1" to "INP". 
c  Subsequent calls to VELOCITY and LCPFCT can work on only portions
c  of the grid, however, to perform restricted integrations on separate
c  line segments.

c  Arguments:
c  (radho) -- I/P real array (INP) : old cell interface positions
c  (radhn) -- I/P real array (INP) : new cell interface positions
c  (i1, inp) -- I/P integer : (first, last) cell interface
c  (alpha) -- I/P integer = (1, 2, 3, 4) => [cartesian, cylindrical, 
c			    spherical, general (user supplied)] geometry
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i1,i1p,i,in,inp,alpha
	parameter (npt=200000)
	real*8 radho(inp),radhn(inp),pi,ftpi

c  /FCT_SCRH/ holds scratch arrays for use by LCPFCT & CNVFCT:
	real*8 scrh(npt),scr1(npt),diff(npt),flxh(npt),fabs(npt),
     &	fsgn(npt),term(npt),terp(npt),lnrhot(npt),lorhot(npt),
     &	rhot(npt),rhotd(npt)
	common /fct_scrh/ scrh,scr1,diff,flxh,fabs,fsgn,term,terp,
     &	lnrhot,lorhot,rhot,rhotd

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo

	data pi,ftpi/0.3141592653589793D+01, 0.4188790204786391D+01/

	i1p = i1 + 1
	in = inp - 1

c  Store the old and new grid interface locations from input and then
c  update the new and average interface and grid coefficients:
	do 1 i = i1,inp
		roh(i) = radho(i)
1		rnh(i) = radhn(i)

c  Select the coordinate system:
	go to (100,200,300,400), alpha

c++++++++++++++++++++++++++++++++++++++
c  Cartesian coordinates:
100	ah(inp) = 1.0d0
	do 101 i = i1,in
		ah(i) = 1.0d0
		lo(i) = roh(i+1) - roh(i)
101		ln(i) = rnh(i+1) - rnh(i)
	go to 500

c++++++++++++++++++++++++++++++++++++++
c  Cylindrical coordinates (radial):
200	diff(i1) = rnh(i1) * rnh(i1)
	scrh(i1) = roh(i1) * roh(i1)
	ah(inp) = pi * (roh(inp) + rnh(inp))
	do 201 i = i1,in
		ah(i) = pi * (roh(i) + rnh(i))
		scrh(i+1) = roh(i+1) * roh(i+1)
		lo(i) = pi * (scrh(i+1) - scrh(i))
		diff(i+1) = rnh(i+1) * rnh(i+1)
201		ln(i) = pi * (diff(i+1) - diff(i))
	go to 500

c++++++++++++++++++++++++++++++++++++++
c  Spherical coordinates (radial):
300	scr1(i1) = roh(i1)**3
	diff(i1) = rnh(i1)**3
	scrh(inp) = (roh(inp) + rnh(inp)) * roh(inp)
	ah(inp) = ftpi * (scrh(inp) + rnh(inp) * rnh(inp))
	do 301 i = i1,in
		scr1(i+1) = roh(i+1)**3
		diff(i+1) = rnh(i+1)**3
		scrh(i) = (roh(i) + rnh(i)) * roh(i)
		ah(i) = ftpi * (scrh(i) + rnh(i) * rnh(i))
		lo(i) = ftpi * (scr1(i+1) - scr1(i))
301		ln(i) = ftpi * (diff(i+1) - diff(i))
	go to 500

c++++++++++++++++++++++++++++++++++++++
c  Special coordinates: Areas and volumes are user-supplied
400	continue

c+++++  Additional system-independent geometric variables  +++++
500	do 501 i = i1,in
	rlo(i) = 1.0d0 / lo(i)
501	rln(i) = 1.0d0 / ln(i)
	lh(i1) = ln(i1)
	rlh(i1) = rln(i1)
	do 502 i = i1p,in
		lh(i) = 0.5d0 * (ln(i) + ln(i-1))
502		rlh(i) = 0.5d0 * (rln(i) + rln(i-1))
	lh(inp) = ln(in)
	rlh(inp) = rln(in)
	do 503 i = i1,inp
503		adugth(i) = ah(i) * (rnh(i) - roh(i))

	return
	END

c********************************************************************

	SUBROUTINE velocity (uh,i1,inp,dt)

c++++++++++++++++++++++++++++++++++++++
c  Description:  Calculates all velocity-dependent coefficients for the
c  LCPFCT & CNVFCT routines. This routine must be called before either
c  LCPFCT or CNVFCT is called. MAKEGRID must be called earlier to set
c  grid and geometry data used here.

c  Arguments:
c  (uh) -- I/P real array (npt) : flow velocity at cell interfaces
c  (i1) -- I/P integer : first cell interface of integration
c  (inp) -- I/P integer : last cell interface = N+1
c  (dt) -- I/P real : stepsize for the time integration
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i1,i1p,i,in,inp
	parameter (npt=200000)
	real*8 uh(inp),dt,rdt,dth,dt2,dt4,one3d,one6th

c  /FCT_SCRH/ holds scratch arrays for use by LCPFCT & CNVFCT:
	real*8 scrh(npt),scr1(npt),diff(npt),flxh(npt),fabs(npt),
     &	fsgn(npt),term(npt),terp(npt),lnrhot(npt),lorhot(npt),
     &	rhot(npt),rhotd(npt)
	common /fct_scrh/ scrh,scr1,diff,flxh,fabs,fsgn,term,terp,
     &	lnrhot,lorhot,rhot,rhotd

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo

c  /FCT_VELO/ holds velocity-dependent flux coefficients:
	real*8 hadudth(npt),nulh(npt),mulh(npt),epsh(npt),vdtodr(npt)
	common /fct_velo/ hadudth,nulh,mulh,epsh,vdtodr

	i1p = i1 + 1
	in = inp - 1

c  Calc. 0.5*interface area * velocity difference * dt (HADUDTH).
c  Next calc. the interface epsilon (EPSH = V*DT/DX). Then find
c  the diffusion (NULH) and antidiffusion (MULH) coefficients. The
c  variation with epsilon gives fourth-order accurate phases when the
c  grid is uniform, the velocity constant, and SCRH is set to zero.
c  With SCRH nonzero (as below) slightly better results are obtained
c  in some of the tests. Optimal performance, of course, depends on
c  the application.

	rdt = 1.0d0 / dt
	dth = 0.5d0 * dt
c???  can the foll. be done in DATA ??
	one6th = 1.0d0 / 6.0d0
	one3d = 1.0d0 / 3.0d0
c?????????????????????????????????????
	do 1 i = i1,inp
		hadudth(i) = dt * ah(i) * uh(i) - adugth(i)
		epsh(i) = hadudth(i) * rlh(i)
		scrh(i) = dmin1(one6th,dabs(epsh(i)))
		scrh(i) = one3d * scrh(i)**2
		hadudth(i) = 0.5d0 * hadudth(i)
		nulh(i) = one6th + one3d * (epsh(i) + scrh(i)) * 
     &			  (epsh(i) - scrh(i))
		mulh(i) = 0.25d0 - 0.5d0 * nulh(i)
		nulh(i) = lh(i) * (nulh(i) + scrh(i))
		mulh(i) = lh(i) * (mulh(i) + scrh(i))
1		diff(i) = uh(i) - rdt * (rnh(i) - roh(i))

c  Now calc. VDTODR for CNVFCT:
	dt2 = 2.0d0 * dt
	dt4 = 4.0d0 * dt
	vdtodr(i1) = dt2 * diff(i1) / (rnh(i1p) - rnh(i1) + 
     &		     roh(i1p) - roh(i1))
	do 2 i = i1p,in
2		vdtodr(i) = dt4 * diff(i) / (rnh(i+1) - rnh(i-1) +
     &					     roh(i+1) - roh(i-1))
	vdtodr(inp) = dt2 * diff(inp) / (rnh(inp) - rnh(in) +
     &					 roh(inp) - roh(in))

	return
	END

c********************************************************************

	SUBROUTINE sources (i1,in,dt,mode,c,d,d1,dn)

c++++++++++++++++++++++++++++++++++++++
c  Description: Accumulates different source terms.

c  Arguments:
c  (i1, in) -- I/P integer : (first, last) cell to be integrated
c  (dt) -- I/P real : stepsize for time integration
c  (mode) -- I/P integer
c		= 1  computes  +DIV(D)
c		= 2  computes  +C*GRAD(D)
c		= 3  adds +D to the sources
c		= 4  +DIV(D) from interface data
c		= 5  +C*GRAD(D) from interface data
c		= 6  +C for list of scalar indices

c  (C) -- I/P real array (NPT) : array of source variables
c  (D) -- I/P real array (NPT) : array of source variables
c  (D1, DN) -- I/P real : (first, last) boundary value of "D".

	implicit real*8 (a-h,o-z)
	integer npt,nindmax,mode,is,i,i1,in,i1p,inp
	parameter (npt=200000, nindmax=150)
	real*8 c(npt),d(npt),dt,dth,dtq,d1,dn

c  /FCT_NDEX/ holds a scalar list of special cell information:
	real*8 scalars(nindmax)
	integer index(nindmax),nind
	common /FCT_NDEX/ nind,index,scalars

c  /FCT_SCRH/ holds scratch arrays for use by LCPFCT & CNVFCT:
	real*8 scrh(npt),scr1(npt),diff(npt),flxh(npt),fabs(npt),
     &	fsgn(npt),term(npt),terp(npt),lnrhot(npt),lorhot(npt),
     &	rhot(npt),rhotd(npt)
	common /fct_scrh/ scrh,scr1,diff,flxh,fabs,fsgn,term,terp,
     &	lnrhot,lorhot,rhot,rhotd

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo

c  /FCT_MISC/ holds the source array and diffusion coefficient:
	real*8 source(npt),diff1
	common /fct_misc/ source,diff1

	i1p = i1 + 1
	inp = in + 1
	dth = 0.5d0 * dt
	dtq = 0.25d0 * dt
	go to (101,202,303,404,505,606), mode

c+++++  [+DIV(D)] is computed conservatively and added to source  +++++
101	scrh(i1) = dt * ah(i1) * d1
	scrh(inp) = dt * ah(inp) * dn
	do 1 i = in,i1p,-1
		scrh(i) = dth * ah(i) * (d(i) + d(i-1))
1		source(i) = source(i) + scrh(i+1) - scrh(i)
	source(i1) = source(i1) + scrh(i1p) - scrh(i1)
	return

c++++  [+C*grad(D)] is computed efficiently and added to SOURCE  +++++
202	scrh(i1) = dth * d1
	scrh(inp) = dth * dn
	do 2 i = in,i1p,-1
		scrh(i) = dtq * (d(i) + d(i-1))
		diff(i) = scrh(i+1) - scrh(i)
2		source(i) = source(i) + c(i) * (ah(i+1)+ah(i)) * diff(i)
	source(i1) = source(i1) + c(i1) * (ah(i1p)+ah(i1)) *
     &		     (scrh(i1p) - scrh(i1))
	return

c+++++  [+D] is added to SOURCE in an explicit formulation  +++++
303	do 3 i = i1,in
3		source(i) = source(i) + dt * lo(i) * d(i)
	return

c+++++  [+div(D)] is computed conservatively from interface data  +++++
404	scrh(inp) = dt * ah(inp) * dn
	scrh(i1) = dt * ah(i1) * d1
	do 4 i = in,i1p,-1
		scrh(i) = dt * ah(i) * d(i)
4		source(i) = source(i) + scrh(i+1) - scrh(i)
	source(i1) = source(i1) + scrh(i1p) - scrh(i1)
	return

c+++++  [+C*grad(D)] is computed using interface data  +++++
505	scrh(i1) = dth * d1
	scrh(inp) = dth * dn
	do 5 i = in,i1p,-1
		scrh(i) = dth * d(i)
		diff(i) = scrh(i+1) - scrh(i)
5		source(i) = source(i) + c(i) * (ah(i+1) + ah(i)) * 
     &			    diff(i)

	source(i1) = source(i1) + c(i1) * (ah(i1p) + ah(i1)) * 
     &		     (scrh(i1p) - scrh(i1))
	return

c+++++  [+C] for source terms only at a list of indices  +++++
606	do 6 is = 1,nind
		i = index(is)
6		source(i) = source(i) + scalars(is)
	return

	END

c********************************************************************

	SUBROUTINE CNVFCT (rhoo,rhon,i1,in,srho1,vrho1,srhon,vrhon,pbc)

c++++++++++++++++++++++++++++++++++++++
c  Description: Solves an advective continuity equation of the form:
c  d(RHO)/dt = -V*grad(RHO) + SOURCES    in the users's choice of 
c  Cartesian, cylindrical or spherical coordinate systems. A facility
c  is included to allow definition of other coordinates. The grid can 
c  be Eulerian, sliding rezone, or Lagrangian, and can be arbitrarily 
c  spaced. The algorithm is a low phase error FCT algorithm, vectorized 
c  and optimized for a combination for speed and flexibility.

c  Arguments: Same as for routine LCPFCT
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i1,in,i1p,inp,i
	real*8 bignum,srho1,vrho1,srhon,vrhon,rho1m,rhonp,rhot1m,rhotnp,
     &	rhotd1m,rhotdnp
	logical pbc
	parameter(npt=200000, bignum=1.0d38)
c  (bignum) = machine-dependent largest number -- set by the user !!

	real*8 rhoo(npt),rhon(npt)

c  /FCT_SCRH/ holds scratch arrays for use by LCPFCT & CNVFCT:
	real*8 scrh(npt),scr1(npt),diff(npt),flxh(npt),fabs(npt),
     &	fsgn(npt),term(npt),terp(npt),lnrhot(npt),lorhot(npt),
     &	rhot(npt),rhotd(npt)
	common /fct_scrh/ scrh,scr1,diff,flxh,fabs,fsgn,term,terp,
     &	lnrhot,lorhot,rhot,rhotd

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo


c  /FCT_VELO/ holds velocity-dependent flux coefficients:
	real*8 hadudth(npt),nulh(npt),mulh(npt),epsh(npt),vdtodr(npt)
	common /fct_velo/ hadudth,nulh,mulh,epsh,vdtodr

c  /FCT_MISC/ holds the source array and diffusion coefficient:
	real*8 source(npt),diff1
	common /fct_misc/ source,diff1

c++++++++++++++++++++++++++++++++++++++
	i1p = i1 + 1
	inp = in + 1

c+++++  Calc. convective and diffusive fluxes  +++++
	if (pbc) then
		rho1m = rhoo(in)
		rhonp = rhoo(i1)
	else
		rho1m = srho1 * rhoo(i1) + vrho1
		rhonp = srhon * rhoo(in) + vrhon
	endif

	diff(i1) = nulh(i1) * (rhoo(i1) - rho1m)
	flxh(i1) = vdtodr(i1) * (rhoo(i1) - rho1m)

	do 1 i = i1p,in
		diff(i) = (rhoo(i) - rhoo(i-1))
		flxh(i) = vdtodr(i) * diff(i)
1		diff(i) = nulh(i) * diff(i)
	diff(inp) = nulh(inp) * (rhonp - rhoo(in))
	flxh(inp) = vdtodr(inp) * (rhonp - rhoo(in))

c  Calc. LORHOT, the transported mass elements, and LNRHOT, the
c  transported and diffused mass elements:
	do 2 i = i1,in
		lorhot(i) = ln(i) * (rhoo(i) - 0.5d0 * 
     &		(flxh(i+1) + flxh(i))) + source(i)

		lnrhot(i) = lorhot(i) + diff(i+1) - diff(i)
		rhot(i) = lorhot(i) * rlo(i)
2		rhotd(i) = lnrhot(i) * rln(i)

c  Evaluate boundary conditions for (rhot, rhotd):
	if (pbc) then
		rhot1m = rhot(in)
		rhotnp = rhot(i1)
		rhotd1m = rhotd(in)
		rhotdnp = rhotd(i1)
	else
		rhot1m = srho1 * rhot(i1) + vrho1
		rhotnp = srhon * rhot(in) + vrhon
		rhotd1m = srho1 * rhotd(i1) + vrho1
		rhotdnp = srhon * rhotd(in) + vrhon
	endif

c  Calc. the transported antidiffusive fluxes and transported and
c  diffused density differences:
	flxh(i1) = mulh(i1) * (rhot(i1) - rhot1m)
	diff(i1) = rhotd(i1) - rhotd1m
	fabs(i1) = dabs(flxh(i1))
	fsgn(i1) = sign(diff1,diff(i1))

	do 3 i = i1p,in
	flxh(i) = mulh(i) * (rhot(i) - rhot(i-1))
3	diff(i) = rhotd(i) - rhotd(i-1)

	flxh(inp) = mulh(inp) * (rhotnp - rhot(in))
	diff(inp) = rhotdnp - rhotd(in)
c++++++++++++++++++++++++++++++++++++++

c  Calc. the magnitude & sign of the antidiffusive flux followed
c  by the flux-limiting changes on the right and left:
	do 4 i = i1,in
	fabs(i+1) = dabs(flxh(i+1))
	fsgn(i+1) = sign(diff1,diff(i+1))
	term(i+1) = fsgn(i+1) * ln(i) * diff(i)
4	terp(i) = fsgn(i) * ln(i) * diff(i+1)

	if (pbc) then
		terp(inp) = terp(i1)
		term(i1) = term(inp)
	else
		terp(inp) = bignum
		term(i1) = bignum
	endif

c  Correct the transported fluxes completely and then calculate the new
c  Flux-Corrected Transport densities:
	flxh(i1) = fsgn(i1) * dmax1(0.0d0,dmin1(term(i1),fabs(i1),
     &						terp(i1)))
	do 5 i = i1,in
	flxh(i+1) = fsgn(i+1) * dmax1(0.0d0,dmin1(term(i+1),fabs(i+1),
     &						  terp(i+1)))
	rhon(i) = rln(i) * (lnrhot(i) + (flxh(i) - flxh(i+1)))
5	source(i) = 0.0d0

	return
	END

c********************************************************************

	SUBROUTINE conserve (rho,i1,in,csum)

c++++++++++++++++++++++++++++++++++++++
c  Description: Computes the ostensibly conserved sum. Beware your
c  boundary conditions and note that only one continuity equation is
c  summed for each call to this routine.

c  Arguments:
c  (rho) -- I/P real array (NPT): cell values for physical var. "rho"
c  (i1, in) -- I/P integer: (first, last) cell to be integrated
c  (csum) -- O/P real: value of the conservation sum of "rho"
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i,i1,in
	parameter (npt=200000)
	real*8 csum,rho(npt)
	
c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo

c  Compute the ostensibly conserved total mas (!!  beware B.C.  !!):
	csum = 0.0d0
	do 80 i = i1,in
80		csum = csum + ln(i) * rho(i)
	return
	END

c********************************************************************

	SUBROUTINE copygrid(mode,i1,in)

c++++++++++++++++++++++++++++++++++++++
c  Description: Makes a complete copy of the grid variables defined by 
c  the most recent call to MAKEGRID from cell I1 to IN, including the 
c  boundary values at interface (IN+1) when the argument MODE=1.
c  When MODE=2, these grid variables are reset from common block 
c  /OLD_GRID/. This routine is used where the same grid is needed 
c  repeatedly after some of the values have been over-written, for 
c  example, by a grid which moves between the halfstep and the whole
c  step.

c  Arguments:
c  (i1, in) -- I/P integer: (first, last) cell index.
c  (mode) -- I/P integer:
c	     = 1  =>  grid variables copied into OLD_GRID
c	     = 2  =>  grid restored from OLD_GRID common.
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i,mode,i1,in
	parameter (npt=200000)

c  /OLD_GRID/ holds geometry, grid, area and volume information:
	real*8 lop(npt),lnp(npt),ahp(npt),rlnp(npt),lhp(npt),rlhp(npt),
     &	rohp(npt),rnhp(npt),adugthp(npt)
	common /old_grid/ lop,lnp,ahp,rlnp,lhp,rlhp,rohp,rnhp,adugthp

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo

	IF (mode.eq.1) THEN
		do 101 i = i1,in
			lop(i) = lo(i)
			lnp(i) = ln(i)
101			rlnp(i) = rln(i)
		do 102 i = i1,in+1
			ahp(i) = ah(i)
			lhp(i) = lh(i)
			rlhp(i) = rlh(i)
			rohp(i) = roh(i)
			rnhp(i) = rnh(i)
102			adugthp(i) = adugth(i)

	ELSE	IF (mode.eq.2) THEN
			do 201 i = i1,in
				lo(i) = lop(i)
				ln(i) = lnp(i)
201				rln(i) = rlnp(i)
			do 202 i = i1,in+1
				ah(i) = ahp(i)
				lh(i) = lhp(i)
				rlh(i) = rlhp(i)
				roh(i) = rohp(i)
				rnh(i) = rnhp(i)
202				adugth(i) = adugthp(i)
		ELSE
			write(*,1001) mode
1001			format(//,' COPYGRID error!!  mode=',i3)
			stop
		ENDIF
	return
	END

c********************************************************************

	BLOCK DATA fctblk

	implicit real*8 (a-h,o-z)
	integer npt
	parameter (npt=200000)
c  /FCT_MISC/ holds the source array and diffusion coefficient:
	real*8 source(npt),diff1
	common /fct_misc/ source,diff1

	data source/npt*0.0d0/,diff1/0.999d0/

	END

c********************************************************************

	SUBROUTINE NEW_GRID (radhn,i1,inp,alpha)

c++++++++++++++++++++++++++++++++++++++
c  Description: Initializes geometry variables and coefficients when 
c  the most recent call to MAKEGRID used the same set of values RADHO 
c  and the new interface locations RADHN are different. NEW_GRID is
c  computationally more efficient than the complete grid procedure in
c  MAKEGRID because several formulae do not need to be recomputed.
c  The grid should generally be defined for the entire number of grid
c  interfaces from 1 to INP. However, subsets of the entire grid may be 
c  reinitialized with care.

c  Arguments:
c  (radhn) -- I/P real array (INP): new cell interface positions.
c  (i1, inp) -- I/P integer (first, last) interface index.
c  (alpha) -- I/P integer = (1, 2, 3, 4) for diff. geometries.
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i1,i1p,i,in,alpha,inp
	parameter (npt=200000)
	real*8 radhn(inp),pi,ftpi

c  /FCT_SCRH/ holds scratch arrays for use by LCPFCT & CNVFCT:
	real*8 scrh(npt),scr1(npt),diff(npt),flxh(npt),fabs(npt),
     &	fsgn(npt),term(npt),terp(npt),lnrhot(npt),lorhot(npt),
     &	rhot(npt),rhotd(npt)
	common /fct_scrh/ scrh,scr1,diff,flxh,fabs,fsgn,term,terp,
     &	lnrhot,lorhot,rhot,rhotd

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo

	data pi,ftpi/0.3141592653589793D+01, 0.4188790204786391D+01/

	i1p = i1 + 1
	in = inp - 1

c  Store the old and new grid interface locations from input and then
c  update the new and average interface and grid coefficients:
	do 1 i = i1,inp
1		rnh(i) = radhn(i)

c  Select the coordinate system:
	go to (100,200,300,400), alpha

c++++++++++++++++++++++++++++++++++++++
c  Cartesian coordinates:
100	ah(inp) = 1.0d0
	do 101 i = i1,in
101		ln(i) = rnh(i+1) - rnh(i)
	go to 500

c++++++++++++++++++++++++++++++++++++++
c  Cylindrical coordinates (radial):
200	diff(i1) = rnh(i1) * rnh(i1)
	ah(inp) = pi * (roh(inp) + rnh(inp))
	do 201 i = i1,in
		ah(i) = pi * (roh(i) + rnh(i))
		diff(i+1) = rnh(i+1) * rnh(i+1)
201		ln(i) = pi * (diff(i+1) - diff(i))
	go to 500

c++++++++++++++++++++++++++++++++++++++
c  Spherical coordinates (radial):
300	diff(i1) = rnh(i1)**3
	scrh(inp) = (roh(inp) + rnh(inp)) * roh(inp)
	ah(inp) = ftpi * (scrh(inp) + rnh(inp) * rnh(inp))
	do 301 i = i1,in
		diff(i+1) = rnh(i+1)**3
		scrh(i) = (roh(i) + rnh(i)) * roh(i)
		ah(i) = ftpi * (scrh(i) + rnh(i) * rnh(i))
301		ln(i) = ftpi * (diff(i+1) - diff(i))
	go to 500

c++++++++++++++++++++++++++++++++++++++
c  Special coordinates: Areas and volumes are user-supplied
400	continue

c+++++  Additional system-independent geometric variables  +++++
500	do 501 i = i1,in
501		rln(i) = 1.0d0 / ln(i)
	lh(i1) = ln(i1)
	rlh(i1) = rln(i1)
	do 502 i = i1p,in
		lh(i) = 0.5d0 * (ln(i) + ln(i-1))
502		rlh(i) = 0.5d0 * (rln(i) + rln(i-1))
	lh(inp) = ln(in)
	rlh(inp) = rln(in)
	do 503 i = i1,inp
503		adugth(i) = ah(i) * (rnh(i) - roh(i))

	return
	END

c********************************************************************

	SUBROUTINE residiff(diffa)

c++++++++++++++++++++++++++++++++++++++
c  Description: Allows the user to give FCT some residual numerical
c  diffusion by making the anti-diffusion coefficient smaller.

c  Arguments:
c  (diffa) -- I/P real: Replacement residual diffusion coefficient.
c			Defaults to 0.999 but could be as high as 1.0000
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt
	real*8 diffa
	parameter (npt=200000)

c  /FCT_MISC/ holds the source array and diffusion coefficient:
	real*8 source(npt),diff1
	common /fct_misc/ source,diff1

	diff1 = diffa

	return
	END

c********************************************************************

	SUBROUTINE set_grid (radr,i1,in)

c++++++++++++++++++++++++++++++++++++++
c  Description: Includes the radial factor in the cell volume for polar
c  coordinates. It must be preceded by a call to MAKE_GRID with ALPHA=1
c  to establish the angular dependence of the cell volumes and areas, 
c  and a call to COPY_GRID to save this angular dependence. The angular
c  coordinate is measured in radians (0 to 2*pi) in cylindrical 
c  coordinates and Cos(theta) (-1 to +1) in spherical coordinates.
c  SET_GRID is called inside the loop over radius in a multi-dimensional
c  model to append the appropriate radial factors when integrating in
c  the angular direction.

c  Arguments:
c  (radr) -- I/P real: radius of cell center.
c  (i1,in) -- I/P integer: (first, last) cell index.
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,i1,i1p,i,in,inp
	real*8 radr
	parameter (npt=200000)

c  /OLD_GRID/ holds geometry, grid, area and volume information:
	real*8 lop(npt),lnp(npt),ahp(npt),rlnp(npt),lhp(npt),rlhp(npt),
     &	rohp(npt),rnhp(npt),adugthp(npt)
	common /old_grid/ lop,lnp,ahp,rlnp,lhp,rlhp,rohp,rnhp,adugthp

c  /FCT_GRID/ holds geometry, grid, area and volume information:
	real*8 lo(npt),ln(npt),ah(npt),rln(npt),lh(npt),rlh(npt),
     &	roh(npt),rnh(npt),adugth(npt),rlo(npt)
	common /fct_grid/ lo,ln,ah,rln,lh,rlh,roh,rnh,adugth,rlo

	i1p = i1 + 1
	inp = in + 1

c  Multiply each volume element by the local radius:
	do 100 i = i1,in
		ln(i) = lnp(i) * radr
100		lo(i) = lop(i) * radr

c  Additional system independent geometric variables:
500	do 501 i = i1,in
501		rln(i) = 1.0d0 / ln(i)
	lh(i1) = ln(i1)
	rlh(i1) = rln(i1)
	do 502 i = i1p,in
		lh(i) = 0.5d0 * (ln(i) + ln(i-1))
502		rlh(i) = 0.5d0 * (rln(i) + rln(i-1))
	lh(inp) = ln(in)
	rlh(inp) = rln(in)

	return
	END

c********************************************************************

	SUBROUTINE zerodiff (ind)

c++++++++++++++++++++++++++++++++++++++
c  Description: Sets the FCT diffusion and anti-diffusion parameters
c  to zero at the specified cell interface to inhibit unwanted
c  diffusion across the interface. This routine is used for inflow
c  and outflow boundary conditions. If argument IND is positive,
c  the coefficients at that particular interface are reset. If IND
c  is negative, the list of NIND indices in INDEX are used to reset
c  that many interface coefficients.

c  Arguments:
c  (ind) -- I/P integer: index of interface to be reset.
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,nindmax,ind,is,i
	parameter (npt=200000, nindmax=150)

c  /FCT_NDEX/ holds a scalar list of special cell information:
	real*8 scalars(nindmax)
	integer index(nindmax),nind
	common /FCT_NDEX/ nind,index,scalars

c  /FCT_VELO/ holds velocity-dependent flux coefficients:
	real*8 hadudth(npt),nulh(npt),mulh(npt),epsh(npt),vdtodr(npt)
	common /fct_velo/ hadudth,nulh,mulh,epsh,vdtodr

	IF (ind.gt.0) THEN
		nulh(ind) = 0.0d0
		mulh(ind) = 0.0d0

c!!  The foll. statement modified !!
c	ELSE IF (ind.le.0) THEN
	ELSE

		if (nind.lt.1 .or. nind.gt.nindmax .or. ind.eq.0) then
		  write(*,*) ' ZERODIFF Error!! IND, NIND = ',ind,nind
		  stop
		endif
		do is = 1,nind
			i = index(is)
			nulh(i) = 0.0d0
			mulh(i) = 0.0d0
		enddo
	ENDIF
	return
	END

c********************************************************************

	SUBROUTINE zeroflux (ind)

c++++++++++++++++++++++++++++++++++++++
c  Description: Sets all velocity dependent FCT parameters to zero at
c  the specified cell interfaces to transport fluxes AND diffusion of 
c  material across the interface. This routine is needed in solid wall
c  boundary conditions. If IND is positive, the coeffs. at that 
c  particular interface are reset. If IND is negative, the list of NIND
c  indices in INDEX are used to reset that many interface coeffs.

c  Arguments:
c  (IND) -- I/P integer: index of interface to be reset.
c++++++++++++++++++++++++++++++++++++++

	implicit real*8 (a-h,o-z)
	integer npt,nindmax,ind,is,i
	parameter (npt=200000, nindmax=150)

c  /FCT_NDEX/ holds a scalar list of special cell information:
	real*8 scalars(nindmax)
	integer index(nindmax),nind
	common /FCT_NDEX/ nind,index,scalars

c  /FCT_VELO/ holds velocity-dependent flux coefficients:
	real*8 hadudth(npt),nulh(npt),mulh(npt),epsh(npt),vdtodr(npt)
	common /fct_velo/ hadudth,nulh,mulh,epsh,vdtodr

	IF (ind.gt.0) THEN
		hadudth(ind) = 0.0d0
		nulh(ind) = 0.0d0
		mulh(ind) = 0.0d0

c!!  The foll. statement modified !!
c	ELSE IF (ind.le.0) THEN
	ELSE
		if (nind.lt.1 .or. nind.gt.nindmax .or. ind.eq.0) then
		  write(*,*) ' ZEROFLUX Error!! IND, NIND = ',ind,nind
		  stop
		endif
		do is = 1,nind
			i = index(is)
			hadudth(i) = 0.0d0
			nulh(i) = 0.0d0
			mulh(i) = 0.0d0
		enddo
	ENDIF
	return
	END
