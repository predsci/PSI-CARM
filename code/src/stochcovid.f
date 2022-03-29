      subroutine stochcovid(nb, Rval, tn,
     $     vecN, vecI0, 
     $     vecInfNess, vecSusTy, vecPS, vecPM, vecPA, matCt, 
     $     matCtClosure, scLim, vecTcalc, vecRtrel, D_E, D_I1, D_I2,
     $      D_HR,D_HD, D_ICU, trig_pres, icu_cap, trickle, isevBchange, 
     $     nAges,nTimes, rtn_inf, rtn_icu, rtn_ded_cum, nReals, iseed,
     $     ndays, yinf, yicu, yded, rt_daily)


      implicit none

      integer nSevClasses, trig_type
      real*8 wl
      parameter (nSevClasses = 4)
      parameter (trig_type = 2)
      parameter (wl = 5.0d0)
      integer nb, nAges, nTimes, isevBchange, nReals, iseed, ndays
      real*8 Rval(nb), tn(nb)
      real*8 save_Rval(nb), save_tn(nb)      
      integer vecN(nAges), vecI0(nAges)
      real*8 vecInfNess(nAges), vecSusTy(nAges)
      real*8 vecPS(nAges), vecPM(nAges), vecPA(nAges)
      real*8 matCt(nAges,nAges),  matCtClosure(nAges,nAges)
      real*8 matCtTmp(nAges,nAges)
      real*8 scLim(2)
      real*8 vecTcalc(nTimes), vecRtrel(nTimes)
      real*8 D_E, D_I1, D_I2, D_HR, D_HD
      real*8 D_ICU(nAges)
      integer trig_pres, trickle
      real*8 icu_cap
      integer totalN, totalI0
      real*8 dt, t_cur, ts
      real*8 vecOne(nAges), vecPF(nAges)
      real*8 vecPMcond(nAges), vecPAcond(nAges)
 
      character sevClassChars(nSevClasses)
      real*8 matPsev(nAges, nSevClasses)
      integer maxICU, capICU
      integer i, j, k, l, ir
      integer ngm_col, ngm_row
      integer rtn_inf(nTimes, nAges, nReals)
      integer rtn_icu(nTimes, nAges, nReals)
      real*8 rtn_ded(nTimes, nAges, nReals)
      real*8 rtn_cded(nTimes, nAges, nReals)       
      real*8 ngm(nAges*nSevClasses,nAges*nSevClasses)
      real*8 EVI(nAges*nSevClasses)
      real*8 EVR(nAges*nSevClasses)
      integer indic(nAges*nSevClasses)
      real*8 VECI(nAges*nSevClasses,nAges*nSevClasses)
      real*8 VECR(nAges*nSevClasses,nAges*nSevClasses)      
      integer nn
      real*8 R0_ngm, beta, beta1, beta2, beta_check, beta_tmp, beta_cur
      
      integer S(nAges), E(nAges), I_1M(nAges), I_2M(nAges)
      integer I_1F(nAges), I_2F(nAges), I_1A(nAges), I_2A(nAges)
      integer I_1S(nAges), I_2S(nAges), R(nAges), ICU(nAges)
      
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

      integer cumICU, sevInICU(nAges), sevOutICU(nAges)
      real*8 vecFoi(nAges), pVecFoi(nAges)

      integer noInf(nAges), noExE(nAges), noEntI1S(nAges)
      integer noEntI1M(nAges), noEntI1A(nAges), noEntI1F(nAges)
      integer noExI_1M(nAges), noExI_2M(nAges)
      integer noExI_1F(nAges), noExI_2F(nAges)
      integer noExI_1S(nAges), noExI_2S(nAges)
      integer noExI_1A(nAges), noExI_2A(nAges)
      integer noExICU(nAges)
      real*8 noDead(nAges)

      integer rtn_inf_cum(nTimes, nAges, nReals)
      real*8 rtn_ded_cum(nTimes, nAges, nReals)
      integer yinf(ndays, nAges, nReals)
      integer yicu(ndays, nAges, nReals)
      integer yded(ndays, nAges, nReals)

      real*8 rt(nTimes)
      real*8 rt_daily(ndays)
      real*8 dinf_int
      
      real*8 u
      integer ignbin
      external ignbin


      call srand(iseed)
      
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
      save_Rval = Rval
      do i = 1, nb
         Rval(i) = Rval(i) / R0_ngm
      enddo

      beta = save_Rval(1) / R0_ngm
      beta_check = save_Rval(1)  / (D_I1 + D_I2)

! convert from absolute days in each R(t) value to relative number of days
      save_tn = tn
      do i = 2, (nb-1)
         tn(i) = tn((i-1)) + tn(i)
      enddo

!     for debugging
      
      beta = Rval(1)
      beta_check = save_Rval(1) / (D_I1 + D_I2)

! Initiate the realisation loop
      
      do ir = 1, nReals
         
        ! Initiate the state variables
        S = vecN - vecI0
        E = 0
        I_1M = vecI0
        I_2M = 0
        I_1F = 0
        I_2F = 0
        I_1A = 0
        I_2A = 0
        I_1S = 0
        I_2S = 0
        R = 0
        ICU = 0

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

        
        rtn_inf(1,:, ir) = vecI0

        rtn_icu(:,:,ir) = 0
        rtn_ded(:,:,ir) = 0
        
        ! Initiate some non-state variables that are needed
        cumICU    = 0
        sevInICU  = 0
        sevOutICU = 0
        
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
     $            dble(trickle)/(dble(totalN) * 7.0 * dt * dble(nAges))


 ! Calculate variable probabilites            
             pVecFoi = 1.0d0 - exp(-dt*vecFoi)

             do i = 1, nAges

                noInf(i) = ignbin(S(i), pVecFoi(i))
               
                noExE(i) = ignbin(E(i), vecProbExE(i))

                noEntI1S(i) = ignbin(noExE(i), vecPS(i))

                noEntI1M(i) = ignbin((noExE(i) - noEntI1S(i)),
     $               vecPMcond(i))
                noEntI1A(i) = ignbin((noExE(i) - noEntI1S(i) -
     $           noEntI1M(i)), vecPAcond(i))
             enddo
             noEntI1F = noExE- noEntI1S- noEntI1M - noEntI1A
             do i = 1, nAges
                noExI_1M(i) = ignbin(I_1M(i), vecProbExI_1M(i))
                noExI_2M(i) = ignbin(I_2M(i), vecProbExI_2M(i))
                noExI_1F(i) = ignbin(I_1F(i), vecProbExI_1F(i))
                noExI_2F(i) = ignbin(I_2F(i), vecProbExI_2F(i))
                noExI_1S(i) = ignbin(I_1S(i), vecProbExI_1S(i))
                noExI_2S(i) = ignbin(I_2S(i), vecProbExI_2S(i))
                noExI_1A(i) = ignbin(I_1A(i), vecProbExI_1A(i))
                noExI_2A(i) = ignbin(I_2A(i), vecProbExI_2A(i))
                noExICU(i)  = ignbin(ICU(i), vecProbExICU(i))
             enddo

! zero before updating properly
             noDead = 0.0d0
             
! Count severe people making it into ICU and those not making
! it into ICU
             capICU = maxICU - sum(ICU)
             if (capICU > sum(ICU)) Then
                sevInICU = sevInICU + noExI_2S
                noDead = dble(noExI_2S * 0.5)
             else if (capICU <= 0) Then 
                sevOutICU = sevOutICU + noExI_2S
                noDead    = noExI_2S * 1.0
             else 
!     TODO age prioritization here
                sevOutICU = sevOutICU + noExI_2S
                noDead    = noExI_2S * 1.0
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

            rtn_inf(j,:, ir) = noInf
            rtn_icu(j,:, ir) = noExI_2S
            rtn_ded(j,:, ir) = noDead
            
            cumICU = cumICU + sum(noExI_2S)

         enddo                  ! End of Loop over Time Steps

      enddo                     ! End of loop over realisations
      
! Make derived outputs
    
        rtn_inf_cum(1,:,:) = rtn_inf(1,:,:)
        rtn_ded_cum(1,:,:) = rtn_ded(1,:,:)
        do i = 2, nTimes
           rtn_inf_cum(i,:,:) = rtn_inf_cum((i-1),:,:) + rtn_inf(i,:,:)
           rtn_ded_cum(i,:,:) = rtn_ded_cum((i-1),:,:) + rtn_ded(i,:,:)
        enddo
        
        call sdaily(nTimes, nAges, nReals, ndays, vecTcalc,
     $       rtn_inf, yinf)

        call sdaily(nTimes, nAges, nReals, ndays, vecTcalc,
     $       rtn_icu, yicu)
       
        call sdaily_dp(nTimes, nAges, nReals, ndays, vecTcalc,
     $       rtn_ded_cum, yded)

        call calc_rt_daily(nTimes, ndays, vecTcalc, rt, rt_daily)
        
      return
      end subroutine stochcovid
