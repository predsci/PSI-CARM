        function calcFit1D(y,gamay,x,wght,ndata)

        implicit none

        integer ndata, i
        real*8 y(ndata),x(ndata), gamay(ndata)
        real*8 wght(ndata)
        real*8 xi, yi,sum,val,CalcFit1D
        real*8 tiny

        tiny = 1.0d0

c x is the simulated data
c y is the base profile

C calculate the P(yi,xi)

        sum = 0.0
        do i=1,ndata

           yi = y(i)
           xi = x(i)+tiny

           val = yi * log(xi) - xi - gamay(i)
           sum = sum + val  * wght(i)

        enddo
        sum = -sum

        CalcFit1D = sum

        return
        end function CalcFit1D
      
      subroutine calc_rt_daily(nTimes, ndays,  xvals, rt, rt_daily)

      implicit none

      integer nTimes, ndays
      real*8 xvals(nTimes)
      real*8 rt(nTimes)
      real*8 rt_daily(ndays)
      real*8 tps((ndays+1))
      real*8 epsilon, t_cur, t_next, tmp
      integer j, i, ind_t, icount 

      epsilon = 1.d-6

      rt_daily = 0.0d0

      tmp = 0.0d0

      do j=1, ndays+1
         tps(j) = dble(j-1)
      enddo

      ind_t = 2
      t_next = tps(ind_t)
      icount = 0
      
      do j = 2, nTimes
         t_cur = xvals(j)
         tmp = tmp + rt(j)
         icount = icount + 1
         if (t_cur > (t_next - epsilon)) Then

            rt_daily(ind_t-1) = tmp/dble(icount)
            tmp = 0.0d0
            icount = 0
            ind_t = ind_t + 1
            ind_t = min(ind_t, (ndays+1))
            t_next = tps(ind_t)
            
         endif
         
      enddo

      return

      end subroutine calc_rt_daily

      
      subroutine ddaily(nTimes, nAges, ndays, xvals, rtn, y)

      implicit none

      integer nTimes, nAges, ndays
      real*8 xvals(nTimes)
      real*8 rtn(nTimes, nAges)
      real*8 y(ndays, nAges)
      real*8 tmp(nAges), tps((ndays+1))
      real*8 epsilon, t_cur, t_next
      integer j, i, ind_t

      epsilon = 1.d-6
      
      y = 0.0d0
      
      tmp = 0.0d0
      t_cur = xvals(1)


      do j=1, ndays+1
         tps(j) = dble(j-1)
      enddo

      ind_t = 2
      t_next = tps(ind_t)
      
      do j = 2, nTimes
         t_cur = xvals(j)
         do i = 1, nAges
            tmp(i)  = tmp(i) + rtn(j, i)
         enddo
        
         if (t_cur > (t_next - epsilon)) Then
            
            do i=1, nAges
               y((ind_t-1),i) = tmp(i)
            enddo
 
            ind_t = ind_t + 1
            ind_t = min(ind_t, (ndays+1))
            t_next = tps(ind_t)
            tmp = 0.0d0
            
            
         end if

 
      enddo
      return
      end subroutine ddaily
      

      
      subroutine sdaily(nTimes, nAges, nReals, ndays, xvals, rtn, y)

      implicit none

      integer nTimes, nAges, nReals, ndays
      real*8 xvals(nTimes)
      integer rtn(nTimes, nAges, nReals)
      integer y(ndays, nAges, nReals)
      integer tmp(nAges)
      real*8 tps((ndays+1))
      real*8  epsilon, t_cur, t_next
      integer k, j, i, ind_t

      epsilon = 1.d-6
      
      y = 0
      
      do j=1, ndays+1
         tps(j) = dble(j-1)
      enddo

      do k = 1, nReals
         tmp = 0
         t_cur = xvals(1)
         
         ind_t = 2
         t_next = tps(ind_t)
         
         do j = 2, nTimes
            t_cur = xvals(j)
            do i = 1, nAges
               tmp(i)  = tmp(i) + rtn(j, i, k)
            enddo
            
            if (t_cur > (t_next - epsilon)) Then
               
               do i=1, nAges
                  y((ind_t-1), i, k) = tmp(i)
               enddo
               
               ind_t = ind_t + 1
               ind_t = min(ind_t, (ndays+1))
               t_next = tps(ind_t)
               tmp = 0
            endif
         enddo
         
      enddo
      
      return
      end subroutine sdaily

      subroutine sdaily_dp(nTimes, nAges, nReals, ndays, xvals, rtn, y)

      implicit none

      integer nTimes, nAges, nReals, ndays
      real*8 xvals(nTimes)
      real*8 rtn(nTimes, nAges, nReals)
      real*8 y(ndays, nAges, nReals)
      real*8 tmp(nAges)
      real*8 tps((ndays+1))
      real*8  epsilon, t_cur, t_next
      integer k, j, i, ind_t

      epsilon = 1.d-6
      
      y = 0
      
      do j=1, ndays+1
         tps(j) = dble(j-1)
      enddo

      do k = 1, nReals
         tmp = 0
         t_cur = xvals(1)
         
         ind_t = 2
         t_next = tps(ind_t)
         
         do j = 2, nTimes
            t_cur = xvals(j)
            do i = 1, nAges
               tmp(i)  = tmp(i) + rtn(j, i, k)
            enddo
            
            if (t_cur > (t_next - epsilon)) Then
               
               do i=1, nAges
                  y((ind_t-1), i, k) = tmp(i)
               enddo
               
               ind_t = ind_t + 1
               ind_t = min(ind_t, (ndays+1))
               t_next = tps(ind_t)
               tmp = 0
            endif
         enddo
         
      enddo
      
      return
      end subroutine sdaily_dp
      
      subroutine ddaily_all(nTimes, nAges, ndays, xvals, rtn, y)

      implicit none

      integer nTimes, nAges, ndays
      real*8 xvals(nTimes)
      real*8 rtn(nTimes, nAges)
      real*8 y(ndays)
      real*8 tmp(nAges), tps((ndays+1))
      real*8 epsilon, t_cur, t_next
      integer j, i, ind_t

      epsilon = 1.d-6
      
      y = 0.0d0
      
      tmp = 0.0d0
      t_cur = xvals(1)


      do j=1, ndays+1
         tps(j) = dble(j-1)
      enddo

      ind_t = 2
      t_next = tps(ind_t)
      
      do j = 2, nTimes
         t_cur = xvals(j)
         do i = 1, nAges
            tmp(i)  = tmp(i) + rtn(j, i)
         enddo
        
         if (t_cur > (t_next - epsilon)) Then
 
            y((ind_t-1)) = sum(tmp)
               
            ind_t = ind_t + 1
            ind_t = min(ind_t, (ndays+1))
            t_next = tps(ind_t)
            tmp = 0.0d0
                       
         end if

      enddo
 
      return
      end subroutine ddaily_all    

      subroutine sdaily_all(nTimes, nAges, ndays, nReals, xvals, rtn, y)

      implicit none

      integer nTimes, nAges, nReals, ndays
      real*8 xvals(nTimes)
      integer rtn(nTimes, nAges, nReals)
      integer y(ndays, nReals)
      integer tmp(nAges)
      real*8 tps((ndays+1))
      real*8 epsilon, t_cur, t_next
      integer k, j, i, ind_t

      epsilon = 1.d-6
      
      y = 0

      do j=1, ndays+1
         tps(j) = dble(j-1)
      enddo      

      do k = 1, nReals 
         tmp = 0
         t_cur = xvals(1)
         ind_t = 2
         t_next = tps(ind_t)
         
         do j = 2, nTimes
            t_cur = xvals(j)
            do i = 1, nAges
               tmp(i)  = tmp(i) + rtn(j, i, k)
            enddo
            
            if (t_cur > (t_next - epsilon)) Then
               
               y((ind_t-1), k) = sum(tmp)
               
               ind_t = ind_t + 1
               ind_t = min(ind_t, (ndays+1))
               t_next = tps(ind_t)
               tmp = 0
               
            end if
         enddo
         
      enddo
      
      return
      end subroutine sdaily_all    
