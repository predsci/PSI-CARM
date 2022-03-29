      subroutine tanhrt(nb, Rval, tn, nTimes, vecTcalc, vecRTrel)
  
  !====================================================================
  ! Hyperbolic Tangent Model for R(t) 
  !
  !====================================================================
      implicit none
      
      integer nb, nTimes
      real*8 Rval(nb), tn(nb)
      real*8 vecTcalc(nTimes), vecRTrel(ntimes)
      real*8 wl
      parameter (wl = 5.0d0)      
      real*8 t_cur, beta_tmp
      integer i, j

 ! convert from absolute days in each R(t) value to relative number of days
 
      do i = 2, (nb-1)
         tn(i) = tn((i-1)) + tn(i)
      enddo
      

      do j = 1, nTimes
         t_cur = vecTcalc(j)
!     Calculate R(t)
         beta_tmp = Rval(1) + Rval(nb)
         do i = 2, nb
            beta_tmp = beta_tmp + (Rval(i) - Rval((i-1))) *
     $           tanh((t_cur-tn((i-1)))/wl)
         enddo
         beta_tmp = beta_tmp * 0.5
         vecRTrel(j) = beta_tmp
      enddo
      
      return
      end subroutine tanhrt
