      subroutine mcmc(nb, Rval, tn, Rval_max, Rval_min, tn_max, tn_min,
     $     Rval_step, tn_step, vecN, vecI0, vecInfNess,
     $     vecSusTy, vecPS, vecPM, vecPA, matCt, matCtClosure,
     $     scLim, vecTcalc, vecRtrel, D_E, D_I1, D_I2, D_HR,D_HD,
     $     D_ICU, trig_pres, icu_cap, trickle, isevBchange, nAges,
     $     nTimes, epi, gamaEpi, wght, ndays, theta, nmcmc, nlines,
     $     iseed)



      implicit none

      integer nSevClasses, trig_type
      parameter (nSevClasses = 4)
      parameter (trig_type = 2)
      integer nb, nAges, nTimes, isevBchange, ndays
      integer nmcmc, nlines, iseed
 
      real*8 Rval(nb), tn(nb)
      real*8 Rval_min(nb), Rval_max(nb)
      real*8 tn_min(nb), tn_max(nb)
      real*8 Rval_step(nb), tn_step(nb)
      real*8 vecN(nAges), vecI0(nAges)
      real*8 vecInfNess(nAges), vecSusTy(nAges)
      real*8 vecPS(nAges), vecPM(nAges), vecPA(nAges)
      real*8 matCt(nAges,nAges),  matCtClosure(nAges,nAges)
      real*8 matCtTmp(nAges,nAges)
      real*8 scLim(2)
      real*8 vecTcalc(nTimes), vecRtrel(nTimes)
      real*8 D_E, D_I1, D_I2, D_HR, D_HD
      real*8 D_ICU(nAges)
      real*8 trig_pres, icu_cap, trickle
      real*8 totalN, totalI0
      real*8 dt, t_cur
      real*8 vecOne(nAges), vecPF(nAges)
      real*8 vecPMcond(nAges), vecPAcond(nAges)
      real*8 Epi(ndays), wght(ndays)
      real*8 gamaEpi(ndays)
      real theta(nlines, (nb*2))
      
      character sevClassChars(nSevClasses)
      real*8 matPsev(nAges, nSevClasses)
      real*8 maxICU, capICU
      integer i, j, k, l
      integer ngm_col, ngm_row
      real*8 rtn_inf(nTimes, nAges)
      real*8 rtn_icu(nTimes, nAges)
      real*8 rtn_ded(nTimes, nAges)
      real*8 rtn_cded(nTimes, nAges)       
      real*8 ngm(nAges*nSevClasses,nAges*nSevClasses)
      real*8 EVI(nAges*nSevClasses)
      real*8 EVR(nAges*nSevClasses)
      integer indic(nAges*nSevClasses)
      real*8 VECI(nAges*nSevClasses,nAges*nSevClasses)
      real*8 VECR(nAges*nSevClasses,nAges*nSevClasses)      
      integer nn
      real*8 R0_ngm, beta, beta1, beta2, beta_check, beta_tmp, beta_cur
      
      real*8 S(nAges), E(nAges), I_1M(nAges), I_2M(nAges)
      real*8 I_1F(nAges), I_2F(nAges), I_1A(nAges), I_2A(nAges)
      real*8 I_1S(nAges), I_2S(nAges), R(nAges), ICU(nAges)
      
      real*8 vecHazExE(nAges), vecHazExICU(nAges)
      real*8 vecHazExI_1M(nAges), vecHazExI_2M(nAges)
      real*8 vecHazExI_1F(nAges), vecHazExI_2F(nAges)
      real*8 vecHazExI_1A(nAges), vecHazExI_2A(nAges)
      real*8 vecHazExI_1S(nAges), vecHazExI_2S(nAges)

      real*8 vecProbExE(nAges),  vecProbExICU(nAges)
      real*8 vecProbExI_1M(nAges), vecProbExI_2M(nAges)
      real*8 vecProbExI_1F(nAges), vecProbExI_2F(nAges)
      real*8 vecProbExI_1S(nAges), vecProbExI_2S(nAges)      
      real*8 vecProbExI_1A(nAges), vecProbExI_2A(nAges)

      real*8 cumICU, sevInICU(nAges), sevOutICU(nAges)
      real*8 vecFoi(nAges), pVecFoi(nAges)

      real*8 noInf(nAges), noExE(nAges), noEntI1S(nAges)
      real*8 noEntI1M(nAges), noEntI1A(nAges), noEntI1F(nAges)
      real*8 noExI_1M(nAges), noExI_2M(nAges)
      real*8 noExI_1F(nAges), noExI_2F(nAges)
      real*8 noExI_1S(nAges), noExI_2S(nAges)
      real*8 noExI_1A(nAges), noExI_2A(nAges)
      real*8 noExICU(nAges), noDead(nAges)

      real*8 rtn_inf_cum(nTimes, nAges)
      real*8 rtn_ded_cum(nTimes, nAges)
      real*8 yinf(ndays, nAges)
      real*8 yicu(ndays, nAges)
      real*8 yded(ndays, nAges)      
      real*8 y(ndays), product, diff(ndays)
      real*8 copy_y(ndays)
      real*8 rt_daily(ndays)      
      real*8 dinf_int

      real*8 cur_rval(nb), cur_tn(nb)
      real*8 updt_rval(nb), updt_tn(nb)
      real*8 save_rval(nb), save_tn(nb)
      real*8 copy_rval(nb), copy_tn(nb)
      real*8 cur_rss, new_rss
      real*8 rv
      real*8 ratio, temp
      
      logical accept
      logical logflag_rval(nb)
      logical logflag_tn(nb)
      integer iaccept, icount, ithin, nbm

      real*8 sum_tn
      real*8 scale, range_max, range_min
      real*8 step_max, step_min, myaccept
      integer ionep, iadapt
      integer ii, jj, kk, ll, mm

      real*8 calcFit1D
      external calaFit1D

      
 ! For an adaptive size MCMC - decide if we need to update step size every 1%
! or 1,000 steps
      scale = 2.0d0
      ionep = int(nmcmc * 0.01)
      ionep = min(1000, ionep)
      range_min = 0.20d0
      range_max = 0.40d0
      step_max = 1.d0
      step_min = 1e-5
      iadapt = 0
      
      nbm = nb - 1
     
      Temp = maxval(epi)/10.0d0
      Temp = max(Temp, 1.0d0)

! Initilize the seed
      
      call srand(iseed)

      logflag_rval(1:nb) = .false.
      logflag_tn(1:nb) = .false.

! ensure that the last tn does not exceed number of days minus what we set to be the minimal value for that tn
      tn(nb) = 0.0
      sum_tn = sum(tn(1:(nb-1)))
      if (sum_tn > (ndays-tn_min(nb))) Then
         tn = tn * dble(ndays-tn_min(nb))/sum(tn(1:(nb-1)))
      endif
      
      call detcovid(nb, Rval, tn,
     $     vecN, vecI0, 
     $     vecInfNess, vecSusTy, vecPS, vecPM, vecPA, matCt, 
     $     matCtClosure, scLim, vecTcalc, vecRtrel, D_E, D_I1, D_I2, 
     $     D_HR,D_HD, D_ICU, trig_pres, icu_cap, trickle, isevBchange, 
     $     nAges, nTimes, rtn_inf, rtn_icu, rtn_ded_cum, ndays,
     $     yinf, yicu, yded, rt_daily)
 
      call ddaily_all(nTimes, nAges, ndays, vecTcalc,
     $     rtn_ded_cum, y)


!     convert back to daily deaths
      
      copy_y = y

      do jj = ndays, 2, -1
         y(jj)   = y(jj) - y(jj-1)
      enddo
      
! calculate goodness of fit
      cur_rss = calcFit1D(epi, gamaEpi, y, wght, ndays)
      
      call dblepr('Initial RSS:', -1, cur_rss, 1)
      cur_rval = Rval
      cur_tn   = tn
      updt_rval = cur_rval
      updt_tn   = cur_tn

      ithin = nmcmc/nlines

      iadapt = 0
      iaccept = 0
      icount = 0
      
      do ii = 1, nmcmc

         save_rval = cur_rval
         save_tn   = cur_tn
     
         call fnProposeParamUpdates(nb, cur_rval,
     $        Rval_min, Rval_max, Rval_step, logflag_rval, updt_rval)
         
         cur_rval =  updt_rval
         
         call fnProposeParamUpdates(nbm,  cur_tn(1:nbm),
     $        tn_min(1:nbm), tn_max(1:nbm), tn_step(1:nbm),
     $        logflag_tn(1:nbm), updt_tn(1:nbm))       

! ensure that the last tn does not exceed number of days
         updt_tn(nb) = 0.0d0

         sum_tn = sum(updt_tn(1:(nb-1)))
         if (sum_tn > (ndays-tn_min(nb))) Then
            updt_tn = updt_tn * dble(ndays-tn_min(nb))
     $           /sum(updt_tn(1:(nb-1)))
         endif
         
c$$$         sum_tn = sum(updt_tn)
c$$$         if (sum_tn > ndays) Then
c$$$            updt_tn = updt_tn * dble(ndays)/sum(updt_tn)
c$$$         endif
         cur_tn = updt_tn
 
         
         call detcovid(nb, cur_rval, cur_tn,
     $     vecN, vecI0, 
     $     vecInfNess, vecSusTy, vecPS, vecPM, vecPA, matCt, 
     $     matCtClosure, scLim, vecTcalc, vecRtrel, D_E, D_I1, D_I2, 
     $     D_HR,D_HD, D_ICU, trig_pres, icu_cap, trickle, isevBchange, 
     $     nAges, nTimes, rtn_inf, rtn_icu, rtn_ded_cum, ndays,
     $     yinf, yicu, yded, rt_daily)
         
         call ddaily_all(nTimes, nAges, ndays, vecTcalc,
     $        rtn_ded_cum, y)

         
!      convert back to daily deaths 

      do jj = ndays, 2, -1
         y(jj)   = y(jj) - y(jj-1)
      enddo

      
!     calculate goodness of fit
      new_rss = calcFit1D(epi, gamaEpi, y, wght, ndays)
      
c$$$       diff = obs - y
c$$$       product = dot_product(diff, diff)
c$$$       new_rss = product/dble(ndays)

       accept = .false.

       rv = rand()

       if (exp(-(new_rss-cur_rss)/Temp)  .gt. rv) Then
          accept = .true.
          iaccept = iaccept + 1
          iadapt = iadapt + 1
          cur_rss = new_rss
          save_rval = cur_rval
          save_tn   = cur_tn
         
       else
          cur_rval = save_rval
          cur_tn   = save_tn
       endif
       if (mod(ii, ithin). eq. 0) Then
          icount = icount + 1
          theta(icount,1:nb)=real(cur_rval(1:nb))
          theta(icount,(nb+1):(2*nb-1))= real(cur_tn(1:(nb-1)))
          theta(icount, (2*nb)) = cur_rss
       endif

c$$$       if (mod(ii,ionep) .eq. 0) Then
c$$$          print*, ii, real(cur_Rval), real(cur_rss)
c$$$       endif
       
        if (mod(ii,ionep) .eq. 0) Then
           myaccept = dble(iadapt)/dble(ionep)
           if (myaccept > range_max) Then
              if((rval_step(1)*scale) < step_max)
     $             rval_step = rval_step * scale
             if((tn_step(1)*scale) < step_max)
     $             tn_step = tn_step * scale               
              endif
              if (myaccept < range_min) Then
             if((rval_step(1)/scale) > step_min)
     $             rval_step= rval_step / scale
             if((tn_step(1)/scale) > step_min)
     $             tn_step = tn_step/scale                    
          endif
          iadapt = 0            !reset acceptance number
c$$$ ! Print information to the screen:
c$$$            call intpr('MCMC Step Number:',-1,ii,1)
c$$$            call dblepr('Theta:',-1, theta(icount,:), (nb*2))
c$$$            call dblepr('accept%',-1,myaccept * 100.0,1)

         endif       
 
      enddo

      call dblepr('Final Accept %:', -1,
     $     dble(iaccept)/dble(nmcmc) * 100.0,1)
            
      return
      end subroutine mcmc
      
c----------------------------------------------------------------
       subroutine fnProposeParamUpdates(nparam,curval,
     $     valmin,valmax,step,logflag,parupdt)

      implicit none
      integer nparam
      real*8 curval(nparam),valmin(nparam),valmax(nparam)
      real*8 parupdt(nparam),step(nparam)
      real*8 x, rv, rtn
      logical logflag(nparam)
      real*8  SR_to_unit, SR_from_unit, ran1
      external SR_to_unit, SR_from_unit, ran1
      integer i


      do i = 1, nparam
 101     continue
          parupdt(i)=curval(i)
          rv = rand()
          rv = (rv - 0.50d0)*step(i)

! convert to a zero - one scale

          x = SR_to_unit(curval(i),valmin(i),valmax(i),
     $         logflag(i))

          x = x + rv

! bring value back to original scale

         rtn = SR_from_unit(x,valmin(i),valmax(i),
     $        logflag(i))

         if (rtn .ge. valmax(i) .or. rtn .le. valmin(i)) go to 101
         parupdt(i) = rtn
         
      enddo
      return
      end subroutine fnProposeParamUpdates
      
c----------------------------------------------------------------

      function SR_to_unit(y,ymin,ymax,logflag)

      implicit none

      real*8 y, ymin,ymax,rtn,SR_to_Unit
      logical logflag


      if (logflag) Then
         rtn = (log10(y) - log10(ymin)) /
     $        (log10(ymax)-log10(ymin))
      else
         rtn = (y - ymin)/(ymax - ymin)
      endif

      SR_to_unit = rtn

      return
      end function SR_to_unit

c----------------------------------------------------------------
      
      function SR_from_unit(x,ymin,ymax,logflag)

      implicit none

      real*8 x,ymin,ymax,rtn,SR_from_unit
      logical logflag


      if (logflag) Then
            rtn = ymin *
     $           10.0**(x*(log10(ymax)-log10(ymin)))
      else
         rtn = ymin + (ymax-ymin)*x
      endif

      SR_from_unit = rtn

      return
      end function SR_from_unit

c----------------------------------------------------------------
      
