      subroutine detcovid(nb, save_Rval, save_tn,
     $     vecN, vecI0, 
     $     vecInfNess, vecSusTy, vecPS, vecPM, vecPA, matCt, 
     $     matCtClosure, scLim, vecTcalc, vecRtrel, D_E, D_I1, D_I2, 
     $     D_HR,D_HD, D_ICU, trig_pres, icu_cap, trickle, isevBchange, 
     $     nAges, nTimes, rtn_inf, rtn_icu, rtn_ded_cum, ndays,
     $     yinf, yicu, yded, rt_daily)



      implicit none

      integer nSevClasses, trig_type
      real*8 wl
      parameter (nSevClasses = 4)
      parameter (trig_type = 2)
      parameter (wl = 5.0d0)
      integer nb, nAges, nTimes, isevBchange, ndays
      real*8 Rval(nb), tn(nb)
      real*8 save_Rval(nb), save_tn(nb)
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
 
      character sevClassChars(nSevClasses)
      real*8 matPsev(nAges, nSevClasses)
      real*8 maxICU, capICU
      integer i, j, k, l
      integer ngm_col, ngm_row
      real*8 rtn_inf(nTimes, nAges)
      real*8 rtn_icu(nTimes, nAges)
      real*8 rtn_ded(nTimes, nAges)
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
      real*8 rt(nTimes)
      real*8 rt_daily(ndays)
      real*8 dinf_int

      
      totalN = sum(vecN)
      totalI0 = sum(vecI0)
      dt = vecTcalc(2) - vecTcalc(1)


      vecOne(1:nAges) = 1.0d0
      vecPF = vecOne - vecPS - vecPM - vecPA

      
      vecPMcond = vecPM / (vecOne - vecPS)
      vecPAcond = vecPA / (vecOne - vecPS - vecPM)


      sevClassChars(1) = "S"
      sevClassChars(2) = "M"
      sevClassChars(3) = "F"
      sevClassChars(4) = "A" 

      matPsev(:,1) = vecPS
      matPsev(:,2) = vecPM
      matPsev(:,3) = vecPF
      matPsev(:,4) = vecPA

      dinf_int = D_I1+D_I2
      if (dt.eq. 1.0d0 .and. D_I1 .eq. 3.0d0 .and. D_I2 .eq. 3.0d0) Then
         dinf_int = 7.060d0
      endif  
! Define Trigger
      maxICU = ceiling(icu_cap * totalN)
        
      do i=1, nSevClasses
         do j =1, nAges
            do k = 1, nSevClasses
               do l = 1, nAges
                  ngm_col = (i-1)*nAges + j
                  ngm_row = (k-1)*nAges + l
                  ngm(ngm_row,ngm_col) = 
     $                 matCt(j,l) * (dinf_int) * matPsev(l,k) *
     $                 vecInfNess(j) * vecSusTy(l)
               enddo
            enddo
         enddo
      enddo

      nn = nAges*nSevClasses
      
      call eigenp(nn, nn, ngm, 24.0d0, EVR, EVI, VECR, VECI, INDIC )
      
      R0_ngm = maxval(evr)
      
      Rval = save_Rval
      tn = save_tn
      
      do i = 1, nb
         Rval(i) = Rval(i) / R0_ngm
      enddo

      beta = save_Rval(1) / R0_ngm
      beta_check = save_Rval(1)  / (D_I1 + D_I2)
! convert from absolute days in each R(t) value to relative number of days
 
      do i = 2, (nb-1)
         tn(i) = tn((i-1)) + tn(i)
      enddo

!     for debugging
      
      beta = Rval(1)
      beta_check = save_Rval(1) / (D_I1 + D_I2)

! Initiate the state variables
        S = vecN - vecI0
        E = 0.0d0
        I_1M = vecI0
        I_2M = 0.0d0
        I_1F = 0.0d0
        I_2F = 0.0d0
        I_1A = 0.0d0
        I_2A = 0.0d0
        I_1S = 0.0d0
        I_2S = 0.0d0
        R = 0.0d0
        ICU = 0.0d0

        ! Initiate constant hazards
        vecHazExE    = 1.0d0/D_E
        vecHazExI_1M = 1.0d0/D_I1
        vecHazExI_2M = 1.0d0/D_I2
        vecHazExI_1F = 1.0d0/D_I1
        vecHazExI_2F = 1.0d0/D_I2
        vecHazExI_1A = 1.0d0/D_I1
        vecHazExI_2A = 1.0d0/D_I2
        vecHazExI_1S = 1.0d0/D_I1
        vecHazExI_2S = 1.0d0/D_I2
        vecHazExICU  = 1.0d0/D_ICU
        
        ! Initiate constant probs
        vecProbExE    = 1.0d0 - exp(-dt * vecHazExE)
        vecProbExI_1M = 1.0d0 - exp(-dt * vecHazExI_1M)
        vecProbExI_2M = 1.0d0 - exp(-dt * vecHazExI_2M)
        vecProbExI_1F = 1.0d0 - exp(-dt * vecHazExI_1F)
        vecProbExI_2F = 1.0d0 - exp(-dt * vecHazExI_2F)
        vecProbExI_1S = 1.0d0 - exp(-dt * vecHazExI_1S)
        vecProbExI_2S = 1.0d0 - exp(-dt * vecHazExI_2S)
        vecProbExI_1A = 1.0d0 - exp(-dt * vecHazExI_1A)
        vecProbExI_2A = 1.0d0 - exp(-dt * vecHazExI_2A)
        vecProbExICU  = 1.0d0 - exp(-dt * vecHazExICU)

        rtn_inf(1,:) = vecI0

        rtn_icu = 0.0d0
        rtn_ded = 0.0d0
        
        ! Initiate some non-state variables that are needed
        cumICU    = 0.0d0
        sevInICU  = 0.0d0 
        sevOutICU = 0.0d0 
        
        do j = 2, nTimes
            t_cur = vecTcalc(j)
!     Calculate R(t)
            beta_tmp = Rval(1) + Rval(nb)
            do i = 2, nb
               beta_tmp = beta_tmp + (Rval(i) - Rval((i-1))) *
     $              tanh((t_cur-tn((i-1)))/wl)
            enddo
            beta_tmp = beta_tmp * 0.5
                             
! Set current beta
            beta_cur = beta_tmp
            
!     Record Rt
            rt(j) = beta_cur * R0_ngm
           
! Set mixing matrix
           
            if ((t_cur < scLim(1)) .or. (t_cur >= scLim(2))) Then
                matCtTmp = matCt
             else 
                matCtTmp = matCtClosure
             endif
             
! Calculate variable hazards

             vecFoi =  beta_cur * vecSusTy * (
     $        (matmul( (I_1M * vecInfNess), matCtTmp))/vecN +
     $            (matmul( (I_2M * vecInfNess), matCtTmp))/vecN +     
     $            (matmul( (I_1F * vecInfNess), matCtTmp))/vecN +
     $            (matmul( (I_2F * vecInfNess), matCtTmp))/vecN +
     $            (matmul( (I_1S * vecInfNess), matCtTmp))/vecN +
     $            (matmul( (I_2S * vecInfNess), matCtTmp))/vecN +              
     $            (matmul( (I_1A * vecInfNess), matCtTmp))/vecN +
     $            (matmul( (I_2A * vecInfNess), matCtTmp))/vecN) +
     $            trickle / (totalN * 7.0d0 * dt * nAges)


 ! Calculate variable probabilites            
             pVecFoi = 1.0d0 - exp(-dt*vecFoi)

             noInf = S * pVecFoi
             noExE = E * vecProbExE
             noEntI1S = noExE * vecPS
             noEntI1M = noExE * vecPM
             noEntI1A = noExE * vecPA
             noEntI1F = noExE * (vecOne - vecPM - vecPS - vecPA)
             noExI_1M = I_1M * vecProbExI_1M
             noExI_2M = I_2M * vecProbExI_2M
             noExI_1F = I_1F * vecProbExI_1F
             noExI_2F = I_2F * vecProbExI_2F
             noExI_1S = I_1S * vecProbExI_1S
             noExI_2S = I_2S * vecProbExI_2S
             noExI_1A = I_1A * vecProbExI_1A
             noExI_2A = I_2A * vecProbExI_2A
             noExICU = ICU * vecProbExICU

! Count severe people making it into ICU and those not making
!     it into ICU
             noDead = 0.0d0
             capICU = maxICU - sum(ICU)
             if (capICU > sum(ICU)) Then
                sevInICU = sevInICU + noExI_2S
                noDead   = noExI_2S * 0.5d0
             else if (capICU <= 0) Then
                sevOutICU = sevOutICU + noExI_2S
                noDead    = noExI_2S * 1.0d0
             else 
!TODO age prioritization here
                sevOutICU = sevOutICU + noExI_2S
                noDead    = noExI_2S * 1.0d0
             endif

            
! Update the state variables
! Problems are probably here
            S = S - noInf
            E = E + noInf - noExE
            I_1M = I_1M + noEntI1M - noExI_1M
            I_2M = I_2M + noExI_1M - noExI_2M
            I_1F = I_1F + noEntI1F - noExI_1F
            I_2F = I_2F + noExI_1F - noExI_2F
            I_1S = I_1S + noEntI1S - noExI_1S
            I_2S = I_2S + noExI_1S - noExI_2S
            I_1A = I_1A + noEntI1A - noExI_1A
            I_2A = I_2A + noExI_1A - noExI_2A
            ICU = ICU + noExI_2S - noExICU
            R = R + noExI_2M + noExI_2F + noExICU + noExI_2A

            rtn_inf(j,:) = noInf
            rtn_icu(j,:) = noExI_2S
            rtn_ded(j,:) = noDead
            
            cumICU = cumICU + sum(noExI_2S)
      
        enddo   ! End of Loop over Time Steps

! Make derived outputs
    
        rtn_inf_cum(1,:) = rtn_inf(1,:)
        rtn_ded_cum(1,:) = rtn_ded(1,:)
        do i = 2, nTimes
           rtn_inf_cum(i,:) = rtn_inf_cum((i-1),:) + rtn_inf(i,:)
           rtn_ded_cum(i,:) = rtn_ded_cum((i-1),:) + rtn_ded(i,:)
        enddo
 
!        call ddaily(nTimes, nAges, ndays, vecTcalc, rtn_ded_cum, y)

        call ddaily(nTimes, nAges, ndays, vecTcalc, rtn_inf,
     $       yinf)

       call ddaily(nTimes, nAges, ndays, vecTcalc, rtn_icu,
     $       yicu)

       call ddaily(nTimes, nAges, ndays, vecTcalc, rtn_ded_cum,
     $       yded)       

       call calc_rt_daily(nTimes, ndays, vecTcalc, rt, rt_daily)
        
      return
      end subroutine detcovid
