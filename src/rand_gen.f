	subroutine GENERATE_RANDOM(NP, RANDOM_NUMBERS)

c*********************************************************************72
c
cc MAIN is the main program for RANDOM_NUMBERS.
c
c  Discussion:
c
c    FORTRAN77 does not have a built in random number generator.
c    One way to deal with this is to write your own.  Here is an example.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2008
c
c  Author:
c
c    John Burkardt
c------------------------------------------------------------------------
      implicit none
      INTEGER NP
      REAL*8 RANDOM_NUMBERS(NP)


      call test01 (NP, RANDOM_NUMBERS )

      end
      

!------------------------------------------------------------

        subroutine GENERATE_RANDOM_2(NP, RANDOM_NUMBERS)

c*********************************************************************72
c
cc MAIN is the main program for RANDOM_NUMBERS.
c
c  Discussion:
c
c    FORTRAN77 does not have a built in random number generator.
c    One way to deal with this is to write your own.  Here is an example.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2008
c
c  Author:
c
c    John Burkardt
c------------------------------------------------------------------------
      implicit none
      INTEGER NP
      REAL*8 RANDOM_NUMBERS(NP)

       call test02 (NP, RANDOM_NUMBERS )

      end



!--------------------------------------------------------
       subroutine test01 (NP, RANDOM_NUMBERS )

      implicit none

      integer i, NP
      REAL*8 RANDOM_NUMBERS(NP)
      double precision r8_uniform_01
      integer seed
      integer seed_init
      double precision u(NP)

      seed_init = 1234567


      seed = seed_init
      
      do i = 1, NP
        u(i) = r8_uniform_01 ( seed )
        RANDOM_NUMBERS(i) =u(i)
      end do

      end




      subroutine test02 (NP, RANDOM_NUMBERS )

      implicit none

      integer NP

      integer i
      double precision r8_uniform_01
      REAL*8  RANDOM_NUMBERS(NP)
      integer seed
      integer seed_in
      integer seed_out
      double precision u(NP)
      double precision u_avg
      double precision u_var

      seed = 123456789

      do i = 1, NP
        seed_in = seed
        u(i) = r8_uniform_01 ( seed )
        seed_out = seed
        RANDOM_NUMBERS(i) =u(i)
      end do

      end






c*********************************************************************
      function r8_uniform_01 ( seed )

c*********************************************************************
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end


