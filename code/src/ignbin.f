      function ignbin( n, pp )
      
c*********************************************************************72
c
cc IGNBIN generates a binomial random deviate.
c
c  Discussion:
c
c    This procedure generates a single random deviate from a binomial
c    distribution whose number of trials is N and whose
c    probability of an event in each trial is P.
c
c    The previous version of this program relied on the assumption that
c    local memory would be preserved between calls.  It set up data
c    one time to be preserved for use over multiple calls.  In the
c    interests of portability, this assumption has been removed, and
c    the "setup" data is recomputed on every call.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Voratas Kachitvichyanukul, Bruce Schmeiser,
c    Binomial Random Variate Generation,
c    Communications of the ACM,
c    Volume 31, Number 2, February 1988, pages 216-222.
c
c  Parameters:
c
c    Input, integer N, the number of binomial trials, from which a
c    random deviate will be generated.
c    0 < N.
cf
c    Input, real*8 PP, the probability of an event in each trial of the
c    binomial distribution from which a random deviate is to be generated.
c    0.0 < PP < 1.0.
c
c    Output, integer IGNBIN, a random deviate from the distribution.
c
      implicit none

      real*8 al
      real*8 alv
      real*8 amaxp
      real*8 c
      real*8 f
      real*8 f1
      real*8 f2
      real*8 ffm
      real*8 fm
      real*8 g
      integer i
      integer ignbin
      integer ix
      integer ix1
      integer k
      integer m
      integer mp
      real*8 pp
      integer n
      real*8 p
      real*8 p1
      real*8 p2
      real*8 p3
      real*8 p4
      real*8 q
      real*8 qn
      real*8 r
      real*8 t
      real*8 u
      real*8 v
      real*8 w
      real*8 w2
      real*8 x
      real*8 x1
      real*8 x2
      real*8 xl
      real*8 xll
      real*8 xlr
      real*8 xm
      real*8 xnp
      real*8 xnpq
      real*8 xr
      real*8 ynorm
      real*8 z
      real*8 z2

      
      if (pp == 0.0D+00 .or. n == 0) then
         ignbin=0
         return
      endif
c$$$      if ( pp < 0.0D+00 .or. 1.0D+00 <= pp ) then
c$$$        write ( *, '(a)' ) ' '
c$$$        write ( *, '(a)' ) 'IGNBIN - Fatal error!'
c$$$        write ( *, '(a)' ) '  PP is out of range.'
c$$$        print*, pp
c$$$        print*,n
c$$$        return
c$$$        !stop 1
c$$$      end if

      p = min ( pp, 1.0D+00 - pp )
      q = 1.0D+00 - p
      xnp = dble ( n ) * p

      if ( xnp .lt. 30.0D+00 ) then
        qn = q ** n
        r = p / q
        g = r * real ( n + 1 )
        go to 20
      end if
c     
c  The calculation of this data was originally intended to be
c  done once, then saved for later calls.  
c
      ffm = xnp + p
      m = ffm
      fm = m
      xnpq = xnp * q
      p1 = nint ( 2.195D+00 * sqrt ( xnpq ) - 4.6D+00 * q ) + 0.5D+00
      xm = fm + 0.5D+00
      xl = xm - p1
      xr = xm + p1
      c = 0.134D+00 + 20.5D+00 / ( 15.3D+00 + fm )
      al = ( ffm - xl ) / ( ffm - xl * p )
      xll = al * ( 1.0D+00 + 0.5D+00 * al )
      al = ( xr - ffm ) / ( xr * q )
      xlr = al * ( 1.0D+00 + 0.5D+00 * al )
      p2 = p1 * ( 1.0D+00 + c + c )
      p3 = p2 + c / xll
      p4 = p3 + c / xlr
c
c  Generate a variate.
c
10    continue

c      u = ran1(iseed) * p4
c      v = ran1(iseed)
      u = rand() * p4
      v = rand()      
c
c  Triangle
c
      if ( u .lt. p1 ) then
        ix = xm - p1 * v + u
        if ( 0.5D+00 .lt. pp ) then
          ix = n - ix
        end if
        ignbin = ix
        return
      end if
c
c  Parallelogram
c
      if ( u .le. p2 ) then

        x = xl + ( u - p1 ) / c
        v = v * c + 1.0D+00 - dabs ( xm - x ) / p1

        if ( v .le. 0.0D+00 .or. 1.0D+00 .lt. v ) then
          go to 10
        end if

        ix = x

      else if ( u .le. p3 ) then

        ix = xl + log ( v ) / xll
        if ( ix .lt. 0 ) then
          go to 10
        end if
        v = v * ( u - p2 ) * xll

      else

        ix = xr - log ( v ) / xlr
        if ( n .lt. ix ) then
          go to 10
        end if
        v = v * ( u - p3 ) * xlr

      end if

      k = abs ( ix - m )

      if ( k .le. 20 .or. xnpq / 2.0 - 1.0 .le. k ) then

        f = 1.0D+00
        r = p / q
        g = ( n + 1 ) * r

        if ( m .lt. ix ) then
          mp = m + 1
          do i = mp, ix
            f = f * ( g / i - r )
          end do
        else if ( ix .lt. m ) then
          ix1 = ix + 1
          do i = ix1, m
            f = f / ( g / i - r )
          end do
        end if

        if ( v .le. f ) then
          if ( 0.5D+00 .lt. pp ) then
            ix = n - ix
          end if
          ignbin = ix
          return
        end if

      else

        amaxp = ( k / xnpq ) * ( ( k * ( k / 3.0D+00
     &    + 0.625D+00 ) + 0.1666666666666D+00 ) / xnpq + 0.5D+00 )
        ynorm = - real ( k * k ) / ( 2.0D+00 * xnpq )
        alv = log ( v )

        if ( alv .lt. ynorm - amaxp ) then
          if ( 0.5D+00 .lt. pp ) then
            ix = n - ix
          end if
          ignbin = ix
        return
        end if

        if ( ynorm + amaxp .lt. alv ) then
          go to 10
        end if

        x1 = real ( ix + 1 )
        f1 = fm + 1.0D+00
        z = real ( n + 1 ) - fm
        w = real ( n - ix + 1 )
        z2 = z * z
        x2 = x1 * x1
        f2 = f1 * f1
        w2 = w * w

        t = xm * log ( f1 / x1 ) + ( n - m + 0.5D+00 ) * log ( z / w ) 
     &    + real ( ix - m ) * log ( w * p / ( x1 * q ))
     &    + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 
     &    - ( 99.0D+00 - 140.0D+00
     &    / f2 ) / f2 ) / f2 ) / f2 ) / f1 / 166320.0D+00
     &    + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 
     &    - ( 99.0D+00 - 140.0D+00
     &    / z2 ) / z2 ) / z2 ) / z2 ) / z / 166320.0D+00
     &    + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 
     &    - ( 99.0D+00 - 140.0D+00
     &    / x2 ) / x2 ) / x2 ) / x2 ) / x1 / 166320.0D+00
     &    + ( 13860.0D+00 - ( 462.0D+00 - ( 132.0D+00 
     &    - ( 99.0D+00 - 140.0D+00 
     &    / w2 ) / w2 ) / w2 ) / w2 ) / w / 166320.0D+00

        if ( alv .le. t ) then
          if ( 0.5D+00 .lt. pp ) then
            ix = n - ix
          end if
          ignbin = ix
          return
        end if

      end if

      go to 10
c
c  Mean less than 30.
c
20    continue

      ix = 0
      f = qn

c     u = ran1(iseed)
      u=rand()

30    continue

      if ( u .lt. f ) then
        if ( 0.5D+00 .lt. pp ) then
          ix = n - ix
        end if
        ignbin = ix
        return
      end if

      if ( ix .le. 110 ) then
        u = u - f
        ix = ix + 1
        f = f * ( g / dble ( ix ) - r )
        go to 30
      end if

      go to 20

      return
      end

      
