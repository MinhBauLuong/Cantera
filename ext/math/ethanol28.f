C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C  A 28-Species Reduced Mechanism for Ethanol/Air with NO Formation
C
C     by Zhaoyu Luo and Tianfeng Lu
C
C       Department of Mechanical Engineering
C       University of Connecticut
C     
C       Email: tlu@engr.uconn.edu
C
C References:
C    "Computational investigations of the effects of thermal 
C     stratification in an ethanol-fuelled HCCI engine," 
C    Fuel, submitted.
C
C    Bhagatwala A., Chen J.H., Lu T.F., 
C    “Direct numerical simulations of HCCI/SACI with ethanol,” 
C    Combust. Flame, 161 (7) 1826-1841, 2014.
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE ETHANOL28  (P, T, Y, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*),WDOT(*)
      DIMENSION RF(180),RB(180),RKLOW(19),C(28)
      DIMENSION XQ(12)
C
      CALL YTCP(P, T, Y, C)
      CALL RATT(T, RF, RB, RKLOW)
      CALL RATX(T, C, RF, RB, RKLOW)
      CALL QSSA(RF, RB, XQ)
      CALL RDOT(RF, RB, WDOT)
C
      END

C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE YTCP (P, T, Y, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), C(*)
      DATA SMALL/1D-50/
C
      C(1) = Y(1)/1.00796998D0
      C(2) = Y(2)/2.01593995D0
      C(3) = Y(3)/1.59994001D1
      C(4) = Y(4)/3.19988003D1
      C(5) = Y(5)/1.70073701D1
      C(6) = Y(6)/1.80153401D1
      C(7) = Y(7)/3.30067703D1
      C(8) = Y(8)/3.40147402D1
      C(9) = Y(9)/2.80105505D1
      C(10) = Y(10)/4.40099506D1
      C(11) = Y(11)/3.00264904D1
      C(12) = Y(12)/6.20252907D1
      C(13) = Y(13)/6.10173208D1
      C(14) = Y(14)/4.80418305D1
      C(15) = Y(15)/4.70338606D1
      C(16) = Y(16)/1.60430303D1
      C(17) = Y(17)/1.50350603D1
      C(18) = Y(18)/2.90621506D1
      C(19) = Y(19)/2.80541806D1
      C(20) = Y(20)/2.70462106D1
      C(21) = Y(21)/4.40535808D1
      C(22) = Y(22)/4.60695207D1
      C(23) = Y(23)/7.7060351D1
      C(24) = Y(24)/1.40066996D1
      C(25) = Y(25)/3.00060997D1
      C(26) = Y(26)/4.60054998D1
      C(27) = Y(27)/4.40127993D1
      C(28) = Y(28)/2.80133991D1
C
      SUM = 0D0
      DO K = 1, 28
         SUM = SUM + C(K)
      ENDDO
      SUM = P/(SUM*T*8.314510D7)
C
      DO K = 1, 28
         C(K) = C(K) * SUM
      ENDDO
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE YTCR (RHO, T, Y, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), C(*)
      DATA SMALL/1D-50/
C
C     h
      C(1) = Y(1)/1.00796998D0
C     h2
      C(2) = Y(2)/2.01593995D0
C     o
      C(3) = Y(3)/1.59994001D1
C     o2
      C(4) = Y(4)/3.19988003D1
C     oh
      C(5) = Y(5)/1.70073701D1
C     h2o
      C(6) = Y(6)/1.80153401D1
C     ho2
      C(7) = Y(7)/3.30067703D1
C     h2o2
      C(8) = Y(8)/3.40147402D1
C     co
      C(9) = Y(9)/2.80105505D1
C     co2
      C(10) = Y(10)/4.40099506D1
C     ch2o
      C(11) = Y(11)/3.00264904D1
C     ho2cho
      C(12) = Y(12)/6.20252907D1
C     o2cho
      C(13) = Y(13)/6.10173208D1
C     ch3o2h
      C(14) = Y(14)/4.80418305D1
C     ch3o2
      C(15) = Y(15)/4.70338606D1
C     ch4
      C(16) = Y(16)/1.60430303D1
C     ch3
      C(17) = Y(17)/1.50350603D1
C     c2h5
      C(18) = Y(18)/2.90621506D1
C     c2h4
      C(19) = Y(19)/2.80541806D1
C     c2h3
      C(20) = Y(20)/2.70462106D1
C     ch3cho
      C(21) = Y(21)/4.40535808D1
C     c2h5oh
      C(22) = Y(22)/4.60695207D1
C     o2c2h4oh
      C(23) = Y(23)/7.7060351D1
C     n
      C(24) = Y(24)/1.40066996D1
C     no
      C(25) = Y(25)/3.00060997D1
C     no2
      C(26) = Y(26)/4.60054998D1
C     n2o
      C(27) = Y(27)/4.40127993D1
C     n2
      C(28) = Y(28)/2.80133991D1
C
      DO K = 1, 28
         C(K) = RHO * C(K)
      ENDDO
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATT (T, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (RU=8.31451D7, SMALL=1.D-200, PATM=1.01325D6)
      DIMENSION RF(*), RB(*), RKLOW(*)
      DIMENSION SMH(40), EG(40)
C
      ALOGT = LOG(T)
      TI = 1D0/T
      TI2 = TI*TI
C
      CALL RDSMH (T, SMH)
      EG(1) = EXP(SMH(1))
      EG(2) = EXP(SMH(2))
      EG(3) = EXP(SMH(3))
      EG(4) = EXP(SMH(4))
      EG(5) = EXP(SMH(5))
      EG(6) = EXP(SMH(6))
      EG(7) = EXP(SMH(7))
      EG(8) = EXP(SMH(8))
      EG(9) = EXP(SMH(9))
      EG(10) = EXP(SMH(10))
      EG(11) = EXP(SMH(11))
      EG(12) = EXP(SMH(12))
      EG(13) = EXP(SMH(13))
      EG(14) = EXP(SMH(14))
      EG(15) = EXP(SMH(15))
      EG(16) = EXP(SMH(16))
      EG(17) = EXP(SMH(17))
      EG(18) = EXP(SMH(18))
      EG(19) = EXP(SMH(19))
      EG(20) = EXP(SMH(20))
      EG(21) = EXP(SMH(21))
      EG(22) = EXP(SMH(22))
      EG(23) = EXP(SMH(23))
      EG(24) = EXP(SMH(24))
      EG(25) = EXP(SMH(25))
      EG(26) = EXP(SMH(26))
      EG(27) = EXP(SMH(27))
      EG(28) = EXP(SMH(28))
      EG(29) = EXP(SMH(29))
      EG(30) = EXP(SMH(30))
      EG(31) = EXP(SMH(31))
      EG(32) = EXP(SMH(32))
      EG(33) = EXP(SMH(33))
      EG(34) = EXP(SMH(34))
      EG(35) = EXP(SMH(35))
      EG(36) = EXP(SMH(36))
      EG(37) = EXP(SMH(37))
      EG(38) = EXP(SMH(38))
      EG(39) = EXP(SMH(39))
      EG(40) = EXP(SMH(40))
      PFAC1 = PATM / (RU*T)
      PFAC2 = PFAC1*PFAC1
      PFAC3 = PFAC2*PFAC1
C
C
      RF(1) = EXP(3.58048786D1 -4.1D-1*ALOGT -8.35339668D3*TI)
      EQK = EG(3)*EG(5)/EG(1)/EG(4)
      RB(1) = RF(1) / MAX(EQK, SMALL)
      RF(2) = EXP(1.08356516D1 +2.67D0*ALOGT -3.16623927D3*TI)
      EQK = EG(1)*EG(5)/EG(2)/EG(3)
      RB(2) = RF(2) / MAX(EQK, SMALL)
      RF(3) = EXP(1.9190789D1 +1.51D0*ALOGT -1.72603317D3*TI)
      EQK = EG(1)*EG(6)/EG(2)/EG(5)
      RB(3) = RF(3) / MAX(EQK, SMALL)
      RF(4) = EXP(1.49040725D1 +2.02D0*ALOGT -6.74310335D3*TI)
      EQK = EG(5)*EG(5)/EG(3)/EG(6)
      RB(4) = RF(4) / MAX(EQK, SMALL)
      RF(5) = EXP(4.52701605D1 -1.4D0*ALOGT -5.25358201D4*TI)
      EQK = EG(1)*EG(1)/EG(2)*PFAC1
      RB(5) = RF(5) / MAX(EQK, SMALL)
      RF(6) = EXP(4.06300863D1 -6.3D-1*ALOGT -5.98324618D4*TI)
      EQK = EG(3)*EG(3)/EG(4)*PFAC1
      RB(6) = RF(6) / MAX(EQK, SMALL)
      RF(7) = EXP(4.14242861D1 -7.4D-1*ALOGT -5.13784218D4*TI)
      EQK = EG(1)*EG(3)/EG(5)*PFAC1
      RB(7) = RF(7) / MAX(EQK, SMALL)
      RF(8) = EXP(5.36049885D1 -1.83D0*ALOGT -5.96311751D4*TI)
      EQK = EG(1)*EG(5)/EG(6)*PFAC1
      RB(8) = RF(8) / MAX(EQK, SMALL)
      RF(9) = EXP(2.80196791D1 +6D-1*ALOGT)
      EQK = EG(7)/EG(1)/EG(4)/PFAC1
      RB(9) = RF(9) / MAX(EQK, SMALL)
      RF(10) = EXP(3.04404238D1 -4.14147317D2*TI)
      EQK = EG(2)*EG(4)/EG(1)/EG(7)
      RB(10) = RF(10) / MAX(EQK, SMALL)
      RF(11) = EXP(3.18907389D1 -1.48448917D2*TI)
      EQK = EG(5)*EG(5)/EG(1)/EG(7)
      RB(11) = RF(11) / MAX(EQK, SMALL)
      RF(12) = 3.25D13
      EQK = EG(4)*EG(5)/EG(3)/EG(7)
      RB(12) = RF(12) / MAX(EQK, SMALL)
      RF(13) = EXP(3.09948627D1 +2.50098684D2*TI)
      EQK = EG(4)*EG(6)/EG(5)/EG(7)
      RB(13) = RF(13) / MAX(EQK, SMALL)
      RF(14) = EXP(3.69688748D1 -3.5D-1*ALOGT -2.50249649D4*TI)
      EQK = EG(7)*EG(7)/EG(4)/EG(8)
      RB(14) = RF(14) / MAX(EQK, SMALL)
      RF(15) = EXP(3.06948792D1 -3.5D-1*ALOGT -1.87599174D4*TI)
      EQK = EG(7)*EG(7)/EG(4)/EG(8)
      RB(15) = RF(15) / MAX(EQK, SMALL)
      RF(16) = EXP(3.33183354D1 -2.43707832D4*TI)
      EQK = EG(5)*EG(5)/EG(8)*PFAC1
      RB(16) = RF(16) / MAX(EQK, SMALL)
      RF(17) = EXP(3.0813233D1 -1.99777017D3*TI)
      EQK = EG(5)*EG(6)/EG(1)/EG(8)
      RB(17) = RF(17) / MAX(EQK, SMALL)
      RF(18) = EXP(2.37913188D1 +1D0*ALOGT -3.01930001D3*TI)
      EQK = EG(2)*EG(7)/EG(1)/EG(8)
      RB(18) = RF(18) / MAX(EQK, SMALL)
      RF(19) = EXP(1.60720517D1 +2D0*ALOGT -1.99777017D3*TI)
      EQK = EG(5)*EG(7)/EG(3)/EG(8)
      RB(19) = RF(19) / MAX(EQK, SMALL)
      RF(20) = EXP(2.83241683D1 -2.1497416D2*TI)
      EQK = EG(6)*EG(7)/EG(5)/EG(8)
      RB(20) = RF(20) / MAX(EQK, SMALL)
      RF(21) = EXP(4.19771599D1 -1.47996022D4*TI)
      EQK = EG(6)*EG(7)/EG(5)/EG(8)
      RB(21) = RF(21) / MAX(EQK, SMALL)
      RF(22) = EXP(2.36136376D1 -1.19966854D3*TI)
      EQK = EG(10)/EG(3)/EG(9)/PFAC1
      RB(22) = RF(22) / MAX(EQK, SMALL)
      RF(23) = EXP(2.76798113D1 -2.1406837D4*TI)
      EQK = EG(3)*EG(10)/EG(4)/EG(9)
      RB(23) = RF(23) / MAX(EQK, SMALL)
      RF(24) = EXP(1.20917835D1 +1.89D0*ALOGT +5.82724901D2*TI)
      EQK = EG(1)*EG(10)/EG(5)/EG(9)
      RB(24) = RF(24) / MAX(EQK, SMALL)
      RF(25) = EXP(3.10355463D1 -1.15739834D4*TI)
      EQK = EG(5)*EG(10)/EG(7)/EG(9)
      RB(25) = RF(25) / MAX(EQK, SMALL)
      RF(26) = EXP(2.68865806D1 +6.6D-1*ALOGT -7.48283185D3*TI)
      EQK = EG(1)*EG(9)/EG(12)*PFAC1
      RB(26) = RF(26) / MAX(EQK, SMALL)
      RF(27) = EXP(2.96565343D1 -2.06318834D2*TI)
      EQK = EG(7)*EG(9)/EG(4)/EG(12)
      RB(27) = RF(27) / MAX(EQK, SMALL)
      RF(28) = 7.34D13
      EQK = EG(2)*EG(9)/EG(1)/EG(12)
      RB(28) = RF(28) / MAX(EQK, SMALL)
      RF(29) = 3.02D13
      EQK = EG(5)*EG(9)/EG(3)/EG(12)
      RB(29) = RF(29) / MAX(EQK, SMALL)
      RF(30) = 3D13
      EQK = EG(1)*EG(10)/EG(3)/EG(12)
      RB(30) = RF(30) / MAX(EQK, SMALL)
      RF(31) = 1.02D14
      EQK = EG(6)*EG(9)/EG(5)/EG(12)
      RB(31) = RF(31) / MAX(EQK, SMALL)
      RF(32) = 2.65D13
      EQK = EG(9)*EG(20)/EG(12)/EG(21)
      RB(32) = RF(32) / MAX(EQK, SMALL)
      RF(33) = EXP(3.3152082D1 -6D-2*ALOGT -7.00477601D3*TI)
      EQK = EG(4)*EG(11)/EG(7)/EG(12)
      RB(33) = RF(33) / MAX(EQK, SMALL)
      RF(34) = 3D13
      EQK = EG(1)*EG(5)*EG(10)/EG(7)/EG(12)*PFAC1
      RB(34) = RF(34) / MAX(EQK, SMALL)
      RF(35) = EXP(6.71085787D1 -4.55D0*ALOGT -2.32989317D4*TI)
      EQK = EG(4)*EG(12)/EG(14)*PFAC1
      RB(35) = RF(35) / MAX(EQK, SMALL)
      RF(36) = EXP(2.83191558D1 -5.86750634D3*TI)
      EQK = EG(12)*EG(13)/EG(11)/EG(14)
      RB(36) = RF(36) / MAX(EQK, SMALL)
      RF(37) = EXP(3.38476272D1 -2.02041492D4*TI)
      EQK = EG(5)*EG(15)/EG(13)*PFAC1
      RB(37) = RF(37) / MAX(EQK, SMALL)
      RF(38) = EXP(3.19485092D1 -1.45932834D4*TI)
      EQK = EG(15)/EG(1)/EG(10)/PFAC1
      RB(38) = RF(38) / MAX(EQK, SMALL)
      RF(39) = EXP(3.21512868D1 +3.7D-1*ALOGT -3.67549454D4*TI)
      EQK = EG(12)*EG(12)/EG(9)/EG(11)
      RB(39) = RF(39) / MAX(EQK, SMALL)
      RF(40) = 3D12
      EQK = EG(2)*EG(9)*EG(9)/EG(12)/EG(12)*PFAC1
      RB(40) = RF(40) / MAX(EQK, SMALL)
      RF(41) = EXP(2.77171988D1 +4.8D-1*ALOGT +1.30836334D2*TI)
      EQK = EG(11)/EG(1)/EG(12)/PFAC1
      RB(41) = RF(41) / MAX(EQK, SMALL)
      RF(42) = EXP(1.75767107D1 +1.5D0*ALOGT -4.00560467D4*TI)
      EQK = EG(11)/EG(2)/EG(9)/PFAC1
      RB(42) = RF(42) / MAX(EQK, SMALL)
      RF(43) = EXP(1.81747802D1 +1.63D0*ALOGT +5.30893584D2*TI)
      EQK = EG(6)*EG(12)/EG(5)/EG(11)
      RB(43) = RF(43) / MAX(EQK, SMALL)
      RF(44) = EXP(1.78655549D1 +1.9D0*ALOGT -1.37881367D3*TI)
      EQK = EG(2)*EG(12)/EG(1)/EG(11)
      RB(44) = RF(44) / MAX(EQK, SMALL)
      RF(45) = EXP(2.2557446D1 +1.15D0*ALOGT -1.13726967D3*TI)
      EQK = EG(5)*EG(12)/EG(3)/EG(11)
      RB(45) = RF(45) / MAX(EQK, SMALL)
      RF(46) = EXP(3.6454499D0 +3.36D0*ALOGT -2.16987027D3*TI)
      EQK = EG(12)*EG(20)/EG(11)/EG(21)
      RB(46) = RF(46) / MAX(EQK, SMALL)
      RF(47) = EXP(-4.94766049D0 +4.52D0*ALOGT -3.31116567D3*TI)
      EQK = EG(8)*EG(12)/EG(7)/EG(11)
      RB(47) = RF(47) / MAX(EQK, SMALL)
      RF(48) = EXP(3.18505288D1 -1.31691802D4*TI)
      EQK = EG(1)*EG(11)/EG(17)*PFAC1
      RB(48) = RF(48) / MAX(EQK, SMALL)
      RF(49) = EXP(-4.2272068D1 +9.5D0*ALOGT +2.76819489D3*TI)
      EQK = EG(7)*EG(11)/EG(4)/EG(17)
      RB(49) = RF(49) / MAX(EQK, SMALL)
      RF(50) = 1.2D13
      EQK = EG(11)*EG(20)/EG(17)/EG(21)
      RB(50) = RF(50) / MAX(EQK, SMALL)
      RF(51) = 2D13
      EQK = EG(2)*EG(11)/EG(1)/EG(17)
      RB(51) = RF(51) / MAX(EQK, SMALL)
      RF(52) = 3.01D11
      EQK = EG(8)*EG(11)/EG(7)/EG(17)
      RB(52) = RF(52) / MAX(EQK, SMALL)
      RF(53) = EXP(2.7014835D1 +4.5D-1*ALOGT -1.81158D3*TI)
      EQK = EG(16)/EG(1)/EG(11)/PFAC1
      RB(53) = RF(53) / MAX(EQK, SMALL)
      RF(54) = EXP(3.4950886D1 -1D0*ALOGT)
      EQK = EG(7)*EG(11)/EG(4)/EG(16)
      RB(54) = RF(54) / MAX(EQK, SMALL)
      RF(55) = EXP(3.3115818D1 -2.52463802D3*TI)
      EQK = EG(7)*EG(11)/EG(4)/EG(16)
      RB(55) = RF(55) / MAX(EQK, SMALL)
      RF(56) = 6D12
      EQK = EG(2)*EG(11)/EG(1)/EG(16)
      RB(56) = RF(56) / MAX(EQK, SMALL)
      RF(57) = 1.2D13
      EQK = EG(8)*EG(11)/EG(7)/EG(16)
      RB(57) = RF(57) / MAX(EQK, SMALL)
      RF(58) = 1.8D14
      EQK = EG(11)*EG(11)/EG(12)/EG(16)
      RB(58) = RF(58) / MAX(EQK, SMALL)
      RF(59) = 2.4D13
      EQK = EG(6)*EG(11)/EG(5)/EG(16)
      RB(59) = RF(59) / MAX(EQK, SMALL)
      RF(60) = 4.2D13
      EQK = EG(5)*EG(11)/EG(3)/EG(16)
      RB(60) = RF(60) / MAX(EQK, SMALL)
      RF(61) = EXP(3.70803784D1 -6D-1*ALOGT -1.92731984D2*TI)
      EQK = EG(20)/EG(1)/EG(21)/PFAC1
      RB(61) = RF(61) / MAX(EQK, SMALL)
      RF(62) = EXP(1.33277502D1 +2.5D0*ALOGT -4.82433819D3*TI)
      EQK = EG(2)*EG(21)/EG(1)/EG(20)
      RB(62) = RF(62) / MAX(EQK, SMALL)
      RF(63) = EXP(1.09733574D1 +2.6D0*ALOGT -1.1020445D3*TI)
      EQK = EG(6)*EG(21)/EG(5)/EG(20)
      RB(63) = RF(63) / MAX(EQK, SMALL)
      RF(64) = EXP(2.07430685D1 +1.5D0*ALOGT -4.32766334D3*TI)
      EQK = EG(5)*EG(21)/EG(3)/EG(20)
      RB(64) = RF(64) / MAX(EQK, SMALL)
      RF(65) = EXP(2.42480273D0 +3.74D0*ALOGT -1.05725822D4*TI)
      EQK = EG(8)*EG(21)/EG(7)/EG(20)
      RB(65) = RF(65) / MAX(EQK, SMALL)
      RF(66) = EXP(1.47156719D1 +2D0*ALOGT -4.16160184D3*TI)
      EQK = EG(21)*EG(21)/EG(20)/EG(22)
      RB(66) = RF(66) / MAX(EQK, SMALL)
      RF(67) = EXP(2.28027074D1 +5D-1*ALOGT +8.83145252D2*TI)
      EQK = EG(2)*EG(11)/EG(5)/EG(21)
      RB(67) = RF(67) / MAX(EQK, SMALL)
      RF(68) = EXP(4.06498002D1 -1.34D0*ALOGT -7.13058018D2*TI)
      EQK = EG(6)*EG(23)/EG(5)/EG(21)
      RB(68) = RF(68) / MAX(EQK, SMALL)
      RF(69) = EXP(1.80558296D1 +1.34D0*ALOGT -5.63602668D3*TI)
      EQK = EG(1)*EG(17)/EG(5)/EG(21)
      RB(69) = RF(69) / MAX(EQK, SMALL)
      RF(70) = EXP(1.72462667D1 +1.6D0*ALOGT -2.2674943D3*TI)
      EQK = EG(1)*EG(16)/EG(5)/EG(21)
      RB(70) = RF(70) / MAX(EQK, SMALL)
      RF(71) = EXP(1.78408622D1 +1.6D0*ALOGT -2.72743434D3*TI)
      EQK = EG(6)*EG(22)/EG(5)/EG(21)
      RB(71) = RF(71) / MAX(EQK, SMALL)
      RF(72) = EXP(2.76310211D1 +2.7D-1*ALOGT +3.45961459D2*TI)
      EQK = EG(5)*EG(17)/EG(7)/EG(21)
      RB(72) = RF(72) / MAX(EQK, SMALL)
      RF(73) = EXP(1.16613455D1 +2.23D0*ALOGT +1.52072077D3*TI)
      EQK = EG(4)*EG(20)/EG(7)/EG(21)
      RB(73) = RF(73) / MAX(EQK, SMALL)
      RF(74) = EXP(3.16456007D1 +5D-2*ALOGT +6.84374668D1*TI)
      EQK = EG(1)*EG(11)/EG(3)/EG(21)
      RB(74) = RF(74) / MAX(EQK, SMALL)
      RF(75) = EXP(2.96520387D1 -1.4251096D4*TI)
      EQK = EG(3)*EG(17)/EG(4)/EG(21)
      RB(75) = RF(75) / MAX(EQK, SMALL)
      RF(76) = EXP(9.71157633D-1 +3.28D0*ALOGT -4.07857109D3*TI)
      EQK = EG(5)*EG(11)/EG(4)/EG(21)
      RB(76) = RF(76) / MAX(EQK, SMALL)
      RF(77) = EXP(2.27789268D1 +9D-1*ALOGT)
      EQK = EG(19)/EG(4)/EG(21)/PFAC1
      RB(77) = RF(77) / MAX(EQK, SMALL)
      RF(78) = RF(36)
      EQK = EG(12)*EG(18)/EG(11)/EG(19)
      RB(78) = RF(78) / MAX(EQK, SMALL)
      RF(79) = EXP(2.59217629D1 -9.29944402D3*TI)
      EQK = EG(18)*EG(21)/EG(19)/EG(20)
      RB(79) = RF(79) / MAX(EQK, SMALL)
      RF(80) = EXP(2.92563324D1 +7.10038718D2*TI)
      EQK = EG(17)*EG(17)/EG(19)/EG(21)
      RB(80) = RF(80) / MAX(EQK, SMALL)
      RF(81) = EXP(2.62326542D1 +7.90050168D2*TI)
      EQK = EG(4)*EG(18)/EG(7)/EG(19)
      RB(81) = RF(81) / MAX(EQK, SMALL)
      RF(82) = EXP(3.71778337D1 -1.61D0*ALOGT -9.35983002D2*TI)
      EQK = EG(4)*EG(17)*EG(17)/EG(19)/EG(19)*PFAC1
      RB(82) = RF(82) / MAX(EQK, SMALL)
      RF(83) = 9.6D13
      EQK = EG(5)*EG(17)/EG(1)/EG(19)
      RB(83) = RF(83) / MAX(EQK, SMALL)
      RF(84) = 3.6D13
      EQK = EG(4)*EG(17)/EG(3)/EG(19)
      RB(84) = RF(84) / MAX(EQK, SMALL)
      RF(85) = EXP(3.4078327D1 -2.1286065D4*TI)
      EQK = EG(5)*EG(17)/EG(18)*PFAC1
      RB(85) = RF(85) / MAX(EQK, SMALL)
      RF(86) = 1D13
      EQK = EG(22)/EG(23)
      RB(86) = RF(86) / MAX(EQK, SMALL)
      RF(87) = EXP(3.04036098D1 +2.86833501D2*TI)
      EQK = EG(21)*EG(21)/EG(20)/EG(23)
      RB(87) = RF(87) / MAX(EQK, SMALL)
      RF(88) = 7D13
      EQK = EG(1)*EG(5)*EG(9)/EG(4)/EG(23)*PFAC1
      RB(88) = RF(88) / MAX(EQK, SMALL)
      RF(89) = 7D13
      EQK = EG(1)*EG(21)/EG(2)/EG(23)
      RB(89) = RF(89) / MAX(EQK, SMALL)
      RF(90) = 3D13
      EQK = EG(22)/EG(23)
      RB(90) = RF(90) / MAX(EQK, SMALL)
      RF(91) = 3D13
      EQK = EG(1)*EG(1)*EG(9)/EG(3)/EG(23)*PFAC1
      RB(91) = RF(91) / MAX(EQK, SMALL)
      RF(92) = 3D13
      EQK = EG(1)*EG(11)/EG(5)/EG(23)
      RB(92) = RF(92) / MAX(EQK, SMALL)
      RF(93) = 3D12
      EQK = EG(9)*EG(11)/EG(10)/EG(23)
      RB(93) = RF(93) / MAX(EQK, SMALL)
      RF(94) = EXP(3.77576522D1 -8D-1*ALOGT)
      EQK = EG(21)/EG(1)/EG(22)/PFAC1
      RB(94) = RF(94) / MAX(EQK, SMALL)
      RF(95) = EXP(2.85064899D1 -7.54825001D2*TI)
      EQK = EG(3)*EG(11)/EG(4)/EG(22)
      RB(95) = RF(95) / MAX(EQK, SMALL)
      RF(96) = 2.41666667D0*RF(95)
      EQK = EG(1)*EG(1)*EG(10)/EG(4)/EG(22)*PFAC1
      RB(96) = RF(96) / MAX(EQK, SMALL)
      RF(97) = 2.08333333D0*RF(95)
      EQK = EG(1)*EG(5)*EG(9)/EG(4)/EG(22)*PFAC1
      RB(97) = RF(97) / MAX(EQK, SMALL)
      RF(98) = 5D13
      EQK = EG(1)*EG(1)*EG(9)/EG(3)/EG(22)*PFAC1
      RB(98) = RF(98) / MAX(EQK, SMALL)
      RF(99) = EXP(2.77089077D1 +4.5D-1*ALOGT -9.16860768D2*TI)
      EQK = EG(24)/EG(1)/EG(25)/PFAC1
      RB(99) = RF(99) / MAX(EQK, SMALL)
      RF(100) = EXP(3.26416564D1 -1.30987299D4*TI)
      EQK = EG(1)*EG(18)/EG(2)/EG(19)
      RB(100) = RF(100) / MAX(EQK, SMALL)
      RF(101) = EXP(2.72539977D1 +1.1D-1*ALOGT +2.16383167D3*TI)
      EQK = EG(25)*EG(25)/EG(24)/EG(26)
      RB(101) = RF(101) / MAX(EQK, SMALL)
      RF(102) = EXP(9.37585481D0 +2.45D0*ALOGT +1.46989589D3*TI)
      EQK = EG(20)*EG(25)/EG(21)/EG(24)
      RB(102) = RF(102) / MAX(EQK, SMALL)
      RF(103) = EXP(3.22047006D1 -1.10707667D2*TI)
      EQK = EG(21)*EG(21)/EG(1)/EG(24)
      RB(103) = RF(103) / MAX(EQK, SMALL)
      RF(104) = 2D12
      EQK = EG(2)*EG(25)/EG(1)/EG(24)
      RB(104) = RF(104) / MAX(EQK, SMALL)
      RF(105) = 1.1D14
      EQK = EG(1)*EG(27)/EG(3)/EG(24)
      RB(105) = RF(105) / MAX(EQK, SMALL)
      RF(106) = 1.1D13
      EQK = EG(5)*EG(31)/EG(7)/EG(24)
      RB(106) = RF(106) / MAX(EQK, SMALL)
      RF(107) = EXP(2.97104627D1 +5.03216668D2*TI)
      EQK = EG(17)*EG(31)/EG(19)/EG(24)
      RB(107) = RF(107) / MAX(EQK, SMALL)
      RF(108) = EXP(2.44798039D1 -5.52028684D2*TI)
      EQK = EG(7)*EG(27)/EG(4)/EG(31)
      RB(108) = RF(108) / MAX(EQK, SMALL)
      RF(109) = EXP(4.63300909D1 -2.02D0*ALOGT -1.04417459D4*TI)
      EQK = EG(11)*EG(21)/EG(31)*PFAC1
      RB(109) = RF(109) / MAX(EQK, SMALL)
      RF(110) = EXP(3.62303471D1 -6.9D-1*ALOGT -1.11865065D4*TI)
      EQK = EG(1)*EG(27)/EG(31)*PFAC1
      RB(110) = RF(110) / MAX(EQK, SMALL)
      RF(111) = EXP(3.35659153D1 -1.01D0*ALOGT -2.38977595D3*TI)
      EQK = EG(7)*EG(25)/EG(4)/EG(24)
      RB(111) = RF(111) / MAX(EQK, SMALL)
      RF(112) = EXP(-9.16290732D-1 +3.88D0*ALOGT -6.85381101D3*TI)
      EQK = EG(7)*EG(25)/EG(4)/EG(24)
      RB(112) = RF(112) / MAX(EQK, SMALL)
      RF(113) = EXP(6.71719992D0 +2.41D0*ALOGT -2.65950009D3*TI)
      EQK = EG(5)*EG(27)/EG(4)/EG(24)
      RB(113) = RF(113) / MAX(EQK, SMALL)
      RF(114) = EXP(3.43762575D1 -7.04503335D3*TI)
      EQK = EG(28)/EG(35)
      RB(114) = RF(114) / MAX(EQK, SMALL)
      RF(115) = 1.17647059D-1*RF(114)
      EQK = EG(29)/EG(35)
      RB(115) = RF(115) / MAX(EQK, SMALL)
      RF(116) = EXP(4.80912325D1 -1.34D0*ALOGT -4.37546893D4*TI)
      EQK = EG(12)*EG(21)/EG(27)*PFAC1
      RB(116) = RF(116) / MAX(EQK, SMALL)
      RF(117) = EXP(3.07964962D1 -1.8327151D3*TI)
      EQK = EG(2)*EG(28)/EG(1)/EG(27)
      RB(117) = RF(117) / MAX(EQK, SMALL)
      RF(118) = EXP(2.94127302D1 -9.40008735D2*TI)
      EQK = EG(5)*EG(28)/EG(3)/EG(27)
      RB(118) = RF(118) / MAX(EQK, SMALL)
      RF(119) = EXP(2.88459339D1 +3.11491117D2*TI)
      EQK = EG(6)*EG(28)/EG(5)/EG(27)
      RB(119) = RF(119) / MAX(EQK, SMALL)
      RF(120) = EXP(3.10355463D1 -1.97009325D4*TI)
      EQK = EG(7)*EG(28)/EG(4)/EG(27)
      RB(120) = RF(120) / MAX(EQK, SMALL)
      RF(121) = EXP(-7.25306646D0 +4.58D0*ALOGT -9.89323969D2*TI)
      EQK = EG(20)*EG(28)/EG(21)/EG(27)
      RB(121) = RF(121) / MAX(EQK, SMALL)
      RF(122) = EXP(2.87329612D1 -5.99834268D3*TI)
      EQK = EG(8)*EG(28)/EG(7)/EG(27)
      RB(122) = RF(122) / MAX(EQK, SMALL)
      RF(123) = RF(122)
      EQK = EG(18)*EG(28)/EG(19)/EG(27)
      RB(123) = RF(123) / MAX(EQK, SMALL)
      RF(124) = EXP(1.20552498D1 +2.4D0*ALOGT -4.10121584D2*TI)
      EQK = EG(6)*EG(29)/EG(5)/EG(27)
      RB(124) = RF(124) / MAX(EQK, SMALL)
      RF(125) = EXP(2.87296334D1 -8.41378268D3*TI)
      EQK = EG(9)*EG(21)/EG(28)*PFAC1
      RB(125) = RF(125) / MAX(EQK, SMALL)
      RF(126) = EXP(3.21252597D1 -6D-1*ALOGT -5.09255268D3*TI)
      EQK = EG(5)*EG(9)*EG(11)/EG(4)/EG(29)*PFAC1
      RB(126) = RF(126) / MAX(EQK, SMALL)
      RF(127) = EXP(2.94360258D1 +2.7D-1*ALOGT -1.40900667D2*TI)
      EQK = EG(25)/EG(1)/EG(26)/PFAC1
      RB(127) = RF(127) / MAX(EQK, SMALL)
      RF(128) = EXP(1.77414365D1 +1.93D0*ALOGT -6.51665585D3*TI)
      EQK = EG(2)*EG(26)/EG(1)/EG(25)
      RB(128) = RF(128) / MAX(EQK, SMALL)
      RF(129) = EXP(1.59630779D1 +1.88D0*ALOGT -9.20886502D1*TI)
      EQK = EG(12)*EG(21)/EG(3)/EG(25)
      RB(129) = RF(129) / MAX(EQK, SMALL)
      RF(130) = 5.82204577D-1*RF(129)
      EQK = EG(1)*EG(29)/EG(3)/EG(25)
      RB(130) = RF(130) / MAX(EQK, SMALL)
      RF(131) = EXP(1.44032972D1 +2D0*ALOGT -1.25804167D3*TI)
      EQK = EG(6)*EG(26)/EG(5)/EG(25)
      RB(131) = RF(131) / MAX(EQK, SMALL)
      RF(132) = EXP(1.89009537D0 +3.7D0*ALOGT -4.78055834D3*TI)
      EQK = EG(20)*EG(26)/EG(21)/EG(25)
      RB(132) = RF(132) / MAX(EQK, SMALL)
      RF(133) = EXP(3.13199006D1 -2.92872101D4*TI)
      EQK = EG(7)*EG(26)/EG(4)/EG(25)
      RB(133) = RF(133) / MAX(EQK, SMALL)
      RF(134) = EXP(2.84330227D1 -8.65029452D3*TI)
      EQK = EG(18)*EG(26)/EG(19)/EG(25)
      RB(134) = RF(134) / MAX(EQK, SMALL)
      RF(135) = 2D13
      EQK = EG(1)*EG(25)/EG(21)/EG(23)
      RB(135) = RF(135) / MAX(EQK, SMALL)
      RF(136) = EXP(6.66124488D1 -5.31D0*ALOGT -3.27090834D3*TI)
      EQK = EG(11)*EG(12)/EG(4)/EG(26)
      RB(136) = RF(136) / MAX(EQK, SMALL)
      RF(137) = EXP(3.39409394D1 -6.1D-1*ALOGT -2.64691967D3*TI)
      EQK = EG(3)*EG(29)/EG(4)/EG(26)
      RB(137) = RF(137) / MAX(EQK, SMALL)
      RF(138) = EXP(5.36526043D1 -1.68D0*ALOGT -4.85100868D4*TI)
      EQK = EG(16)*EG(21)/EG(30)*PFAC1
      RB(138) = RF(138) / MAX(EQK, SMALL)
      RF(139) = EXP(5.38349259D1 -1.62D0*ALOGT -5.00901871D4*TI)
      EQK = EG(5)*EG(24)/EG(30)*PFAC1
      RB(139) = RF(139) / MAX(EQK, SMALL)
      RF(140) = EXP(1.17905572D1 +2.52D0*ALOGT -3.05251231D4*TI)
      EQK = EG(6)*EG(25)/EG(30)*PFAC1
      RB(140) = RF(140) / MAX(EQK, SMALL)
      RF(141) = EXP(2.73080572D1 +1D-1*ALOGT -4.57977489D4*TI)
      EQK = EG(2)*EG(27)/EG(30)*PFAC1
      RB(141) = RF(141) / MAX(EQK, SMALL)
      RF(142) = EXP(3.06267534D1 -2.65698401D4*TI)
      EQK = EG(7)*EG(32)/EG(4)/EG(30)
      RB(142) = RF(142) / MAX(EQK, SMALL)
      RF(143) = EXP(3.03390713D1 -2.52363159D4*TI)
      EQK = EG(7)*EG(33)/EG(4)/EG(30)
      RB(143) = RF(143) / MAX(EQK, SMALL)
      RF(144) = EXP(2.59217629D1 +4D-1*ALOGT -3.60806351D2*TI)
      EQK = EG(6)*EG(32)/EG(5)/EG(30)
      RB(144) = RF(144) / MAX(EQK, SMALL)
      RF(145) = EXP(2.4741449D1 +5D-1*ALOGT +1.91222334D2*TI)
      EQK = EG(6)*EG(33)/EG(5)/EG(30)
      RB(145) = RF(145) / MAX(EQK, SMALL)
      RF(146) = EXP(2.3431316D1 +8D-1*ALOGT -1.27515104D3*TI)
      EQK = EG(6)*EG(31)/EG(5)/EG(30)
      RB(146) = RF(146) / MAX(EQK, SMALL)
      RF(147) = EXP(7.53902706D0 +3.2D0*ALOGT -3.59799917D3*TI)
      EQK = EG(2)*EG(32)/EG(1)/EG(30)
      RB(147) = RF(147) / MAX(EQK, SMALL)
      RF(148) = EXP(1.20951411D1 +2.53D0*ALOGT -1.721001D3*TI)
      EQK = EG(2)*EG(33)/EG(1)/EG(30)
      RB(148) = RF(148) / MAX(EQK, SMALL)
      RF(149) = EXP(1.08893043D1 +2.53D0*ALOGT -2.21666942D3*TI)
      EQK = EG(2)*EG(31)/EG(1)/EG(30)
      RB(149) = RF(149) / MAX(EQK, SMALL)
      RF(150) = EXP(1.00770206D1 +2.55D0*ALOGT -8.29804285D3*TI)
      EQK = EG(8)*EG(32)/EG(7)/EG(30)
      RB(150) = RF(150) / MAX(EQK, SMALL)
      RF(151) = EXP(2.94227806D1 -8.05146668D3*TI)
      EQK = EG(8)*EG(33)/EG(7)/EG(30)
      RB(151) = RF(151) / MAX(EQK, SMALL)
      RF(152) = EXP(2.85473118D1 -1.20772D4*TI)
      EQK = EG(8)*EG(31)/EG(7)/EG(30)
      RB(152) = RF(152) / MAX(EQK, SMALL)
      RF(153) = EXP(9.41735454D0 +2.55D0*ALOGT -7.92566252D3*TI)
      EQK = EG(18)*EG(32)/EG(19)/EG(30)
      RB(153) = RF(153) / MAX(EQK, SMALL)
      RF(154) = EXP(9.01188943D0 +2.55D0*ALOGT -5.40957918D3*TI)
      EQK = EG(18)*EG(33)/EG(19)/EG(30)
      RB(154) = RF(154) / MAX(EQK, SMALL)
      RF(155) = RF(152)
      EQK = EG(18)*EG(31)/EG(19)/EG(30)
      RB(155) = RF(155) / MAX(EQK, SMALL)
      RF(156) = EXP(6.87626461D0 +3.23D0*ALOGT -2.34398324D3*TI)
      EQK = EG(5)*EG(32)/EG(3)/EG(30)
      RB(156) = RF(156) / MAX(EQK, SMALL)
      RF(157) = EXP(1.1884489D1 +2.47D0*ALOGT -4.40817801D2*TI)
      EQK = EG(5)*EG(33)/EG(3)/EG(30)
      RB(157) = RF(157) / MAX(EQK, SMALL)
      RF(158) = EXP(-6.52931884D0 +4.73D0*ALOGT -8.69055185D2*TI)
      EQK = EG(5)*EG(31)/EG(3)/EG(30)
      RB(158) = RF(158) / MAX(EQK, SMALL)
      RF(159) = EXP(5.79909265D0 +3.3D0*ALOGT -6.18453285D3*TI)
      EQK = EG(20)*EG(32)/EG(21)/EG(30)
      RB(159) = RF(159) / MAX(EQK, SMALL)
      RF(160) = EXP(2.99222613D0 +3.37D0*ALOGT -3.84155604D3*TI)
      EQK = EG(20)*EG(33)/EG(21)/EG(30)
      RB(160) = RF(160) / MAX(EQK, SMALL)
      RF(161) = EXP(7.10495819D-1 +3.57D0*ALOGT -3.88533589D3*TI)
      EQK = EG(20)*EG(31)/EG(21)/EG(30)
      RB(161) = RF(161) / MAX(EQK, SMALL)
      RF(162) = EXP(5.76105563D1 -3.99D0*ALOGT -1.52927545D4*TI)
      EQK = EG(5)*EG(25)/EG(32)*PFAC1
      RB(162) = RF(162) / MAX(EQK, SMALL)
      RF(163) = EXP(3.22361913D1 -1.25804167D4*TI)
      EQK = EG(1)*EG(27)/EG(33)*PFAC1
      RB(163) = RF(163) / MAX(EQK, SMALL)
      RF(164) = EXP(3.8202338D1 -1D0*ALOGT -1.50965D4*TI)
      EQK = EG(4)*EG(32)/EG(34)*PFAC1
      RB(164) = RF(164) / MAX(EQK, SMALL)
      RF(165) = EXP(2.18627001D1 -9.51079502D3*TI)
      EQK = EG(5)*EG(11)*EG(11)/EG(34)*PFAC2
      RB(165) = RF(165) / MAX(EQK, SMALL)
      RF(166) = EXP(1.51531397D1 +2D0*ALOGT -8.25778552D2*TI)
      EQK = EG(7)*EG(27)/EG(4)/EG(33)
      RB(166) = RF(166) / MAX(EQK, SMALL)
      RF(167) = 3.8D13
      EQK = EG(1)*EG(37)/EG(5)/EG(36)
      RB(167) = RF(167) / MAX(EQK, SMALL)
      RF(168) = EXP(2.25795638D1 +1D0*ALOGT -3.16020067D3*TI)
      EQK = EG(3)*EG(37)/EG(4)/EG(36)
      RB(168) = RF(168) / MAX(EQK, SMALL)
      RF(169) = 2.1D13
      EQK = EG(3)*EG(40)/EG(36)/EG(37)
      RB(169) = RF(169) / MAX(EQK, SMALL)
      RF(170) = EXP(3.48011407D1 -7.5D-1*ALOGT)
      EQK = EG(38)/EG(3)/EG(37)/PFAC1
      RB(170) = RF(170) / MAX(EQK, SMALL)
      RF(171) = 7.26643599D-2*RF(13)
      EQK = EG(5)*EG(38)/EG(7)/EG(37)
      RB(171) = RF(171) / MAX(EQK, SMALL)
      RF(172) = EXP(3.24985556D1 -1.82164434D2*TI)
      EQK = EG(5)*EG(37)/EG(1)/EG(38)
      RB(172) = RF(172) / MAX(EQK, SMALL)
      RF(173) = EXP(3.23315015D1 -5.2D-1*ALOGT)
      EQK = EG(4)*EG(37)/EG(3)/EG(38)
      RB(173) = RF(173) / MAX(EQK, SMALL)
      RF(174) = EXP(2.91350985D1 -1.38882768D4*TI)
      EQK = EG(4)*EG(37)*EG(37)/EG(38)/EG(38)*PFAC1
      RB(174) = RF(174) / MAX(EQK, SMALL)
      RF(175) = EXP(2.78933854D1 -3.14862669D4*TI)
      EQK = EG(3)*EG(40)/EG(39)*PFAC1
      RB(175) = RF(175) / MAX(EQK, SMALL)
      RF(176) = EXP(1.79743936D1 +1.835D0*ALOGT -6.78939928D3*TI)
      EQK = EG(5)*EG(40)/EG(1)/EG(39)
      RB(176) = RF(176) / MAX(EQK, SMALL)
      RF(177) = EXP(3.21528097D1 -1.39285341D4*TI)
      EQK = EG(37)*EG(37)/EG(3)/EG(39)
      RB(177) = RF(177) / MAX(EQK, SMALL)
      RF(178) = EXP(2.89393539D1 -8.01926082D3*TI)
      EQK = EG(4)*EG(40)/EG(3)/EG(39)
      RB(178) = RF(178) / MAX(EQK, SMALL)
      RF(179) = EXP(-4.34280592D0 +4.72D0*ALOGT -1.83976014D4*TI)
      EQK = EG(7)*EG(40)/EG(5)/EG(39)
      RB(179) = RF(179) / MAX(EQK, SMALL)
      RF(180) = EXP(1.31806323D1 +2.23D0*ALOGT -2.32888674D4*TI)
      EQK = EG(38)*EG(40)/EG(37)/EG(39)
      RB(180) = RF(180) / MAX(EQK, SMALL)
C
      RKLOW(1) = EXP(3.80889683D1 -4.11D-1*ALOGT +5.61086584D2*TI)
      RKLOW(2) = EXP(3.93279334D1 -2.28963584D4*TI)
      RKLOW(3) = EXP(5.55621468D1 -2.788D0*ALOGT -2.10898105D3*TI)
      RKLOW(4) = EXP(5.55621468D1 -2.57D0*ALOGT -7.17083751D2*TI)
      RKLOW(5) = EXP(6.37931383D1 -3.42D0*ALOGT -4.24453195D4*TI)
      RKLOW(6) = EXP(5.81889602D1 -3D0*ALOGT -1.22316875D4*TI)
      RKLOW(7) = EXP(7.39217399D1 -4.82D0*ALOGT -3.28600484D3*TI)
      RKLOW(8) = EXP(7.66692127D1 -4.76D0*ALOGT -1.22986154D3*TI)
      RKLOW(9) = EXP(5.71862909D1 -3D0*ALOGT)
      RKLOW(10) = EXP(6.33329483D1 -3.14D0*ALOGT -6.18956501D2*TI)
      RKLOW(11) = EXP(9.68908955D1 -7.62D0*ALOGT -3.50742017D3*TI)
      RKLOW(12) = EXP(3.4721098D1 -6.29926625D3*TI)
      RKLOW(13) = EXP(6.9414025D1 -3.86D0*ALOGT -1.67067934D3*TI)
      RKLOW(14) = EXP(1.96854356D2 -1.884D1*ALOGT -5.69138051D4*TI)
      RKLOW(15) = EXP(1.97350932D2 -1.88D1*ALOGT -5.97670436D4*TI)
      RKLOW(16) = EXP(1.27770351D2 -1.092D1*ALOGT -3.15235049D4*TI)
      RKLOW(17) = EXP(2.01820052D2 -1.942D1*ALOGT -5.81617824D4*TI)
      RKLOW(18) = EXP(5.6813851D1 -2.87D0*ALOGT -7.79985835D2*TI)
      RKLOW(19) = EXP(3.36224857D1 -2.84820634D4*TI)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDSMH  (T, SMH)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION SMH(*), TN(5)
C
      TLOG = LOG(T)
      TI = 1D0/T
C
      TN(1) = TLOG - 1D0
      TN(2) = T
      TN(3) = TN(2)*T
      TN(4) = TN(3)*T
      TN(5) = TN(4)*T
C h
      IF (T .GT. 1000) THEN
      SMH(1) = -4.601176D-1 -2.547163D4*TI 
     *         +2.5D0*TN(1) 
      ELSE
      SMH(1) = -4.601176D-1 -2.547163D4*TI 
     *         +2.5D0*TN(1) 
      ENDIF
C h2
      IF (T .GT. 1000) THEN
      SMH(2) = -1.35511D0 +8.35034D2*TI 
     *         +2.991423D0*TN(1) +3.500322D-4*TN(2) 
     *         -9.389715D-9*TN(3) -7.69298167D-13*TN(4) 
     *         +7.91376D-17*TN(5) 
      ELSE
      SMH(2) = -3.294094D0 +1.012521D3*TI 
     *         +3.298124D0*TN(1) +4.124721D-4*TN(2) 
     *         -1.35716917D-7*TN(3) -7.896195D-12*TN(4) 
     *         +2.067436D-14*TN(5) 
      ENDIF
C o
      IF (T .GT. 1000) THEN
      SMH(3) = 4.920308D0 -2.92308D4*TI 
     *         +2.54206D0*TN(1) -1.377531D-5*TN(2) 
     *         -5.17133833D-10*TN(3) +3.79255583D-13*TN(4) 
     *         -2.184026D-17*TN(5) 
      ELSE
      SMH(3) = 2.963995D0 -2.914764D4*TI 
     *         +2.946429D0*TN(1) -8.19083D-4*TN(2) 
     *         +4.03505333D-7*TN(3) -1.3357025D-10*TN(4) 
     *         +1.945348D-14*TN(5) 
      ENDIF
C o2
      IF (T .GT. 1000) THEN
      SMH(4) = 3.189166D0 +1.23393D3*TI 
     *         +3.697578D0*TN(1) +3.0675985D-4*TN(2) 
     *         -2.09807D-8*TN(3) +1.47940083D-12*TN(4) 
     *         -5.682175D-17*TN(5) 
      ELSE
      SMH(4) = 6.034738D0 +1.005249D3*TI 
     *         +3.212936D0*TN(1) +5.63743D-4*TN(2) 
     *         -9.59358333D-8*TN(3) +1.0948975D-10*TN(4) 
     *         -4.384277D-14*TN(5) 
      ENDIF
C oh
      IF (T .GT. 1710) THEN
      SMH(5) = 5.78756825D0 -3.6994972D3*TI 
     *         +2.8537604D0*TN(1) +5.1497167D-4*TN(2) 
     *         -3.87777462D-8*TN(3) +1.6145892D-12*TN(4) 
     *         -1.57879924D-17*TN(5) 
      ELSE
      SMH(5) = 2.54433372D0 -3.45264448D3*TI 
     *         +3.41896226D0*TN(1) +1.596279D-4*TN(2) 
     *         -5.13821195D-8*TN(3) +3.03672912D-11*TN(4) 
     *         -5.00977395D-15*TN(5) 
      ENDIF
C h2o
      IF (T .GT. 1000) THEN
      SMH(6) = 6.862817D0 +2.989921D4*TI 
     *         +2.672146D0*TN(1) +1.5281465D-3*TN(2) 
     *         -1.45504333D-7*TN(3) +1.00083D-11*TN(4) 
     *         -3.195809D-16*TN(5) 
      ELSE
      SMH(6) = 2.590233D0 +3.020811D4*TI 
     *         +3.386842D0*TN(1) +1.737491D-3*TN(2) 
     *         -1.059116D-6*TN(3) +5.80715083D-10*TN(4) 
     *         -1.253294D-13*TN(5) 
      ENDIF
C ho2
      IF (T .GT. 1000) THEN
      SMH(7) = 3.78510215D0 -1.11856713D2*TI 
     *         +4.0172109D0*TN(1) +1.11991007D-3*TN(2) 
     *         -1.05609692D-7*TN(3) +9.52053083D-12*TN(4) 
     *         -5.39542675D-16*TN(5) 
      ELSE
      SMH(7) = 3.71666245D0 -2.9480804D2*TI 
     *         +4.30179801D0*TN(1) -2.37456026D-3*TN(2) 
     *         +3.52638152D-6*TN(3) -2.02303245D-9*TN(4) 
     *         +4.64612562D-13*TN(5) 
      ENDIF
C h2o2
      IF (T .GT. 1000) THEN
      SMH(8) = 5.01137D-1 +1.800696D4*TI 
     *         +4.573167D0*TN(1) +2.168068D-3*TN(2) 
     *         -2.457815D-7*TN(3) +1.95742D-11*TN(4) 
     *         -7.15827D-16*TN(5) 
      ELSE
      SMH(8) = 6.785363D0 +1.766315D4*TI 
     *         +3.388754D0*TN(1) +3.284613D-3*TN(2) 
     *         -2.47502167D-8*TN(3) -3.85483833D-10*TN(4) 
     *         +1.2357575D-13*TN(5) 
      ENDIF
C co
      IF (T .GT. 1429) THEN
      SMH(9) = 5.71725177D0 +1.42718539D4*TI 
     *         +3.1121689D0*TN(1) +5.79741415D-4*TN(2) 
     *         -5.64133937D-8*TN(3) +3.67835915D-12*TN(4) 
     *         -1.06431114D-16*TN(5) 
      ELSE
      SMH(9) = 5.33277914D0 +1.42869054D4*TI 
     *         +3.19036352D0*TN(1) +4.47209986D-4*TN(2) 
     *         -5.41545938D-9*TN(3) -8.71666392D-12*TN(4) 
     *         +1.20982847D-15*TN(5) 
      ENDIF
C co2
      IF (T .GT. 1380) THEN
      SMH(10) = -5.18289303D0 +4.93178953D4*TI 
     *         +5.18953018D0*TN(1) +1.03003238D-3*TN(2) 
     *         -1.22262554D-7*TN(3) +9.7503645D-12*TN(4) 
     *         -3.45864607D-16*TN(5) 
      ELSE
      SMH(10) = 8.81141041D0 +4.8416283D4*TI 
     *         +2.5793049D0*TN(1) +4.12342494D-3*TN(2) 
     *         -1.07119341D-6*TN(3) +2.1219752D-10*TN(4) 
     *         -2.06015222D-14*TN(5) 
      ENDIF
C ch2o
      IF (T .GT. 1486) THEN
      SMH(11) = 1.06525547D0 +1.49287258D4*TI 
     *         +4.02068394D0*TN(1) +2.54951708D-3*TN(2) 
     *         -2.940508D-7*TN(3) +2.30021566D-11*TN(4) 
     *         -8.0499021D-16*TN(5) 
      ELSE
      SMH(11) = 8.10120233D0 +1.41188397D4*TI 
     *         +3.00754197D0*TN(1) +1.52364748D-3*TN(2) 
     *         +8.75182077D-7*TN(3) -4.26682734D-10*TN(4) 
     *         +6.35668975D-14*TN(5) 
      ENDIF
C hco
      IF (T .GT. 1690) THEN
      SMH(12) = 6.24593456D0 -3.97409684D3*TI 
     *         +3.44148164D0*TN(1) +1.76078859D-3*TN(2) 
     *         -2.0689353D-7*TN(3) +1.64440537D-11*TN(4) 
     *         -5.8269308D-16*TN(5) 
      ELSE
      SMH(12) = 4.94843165D0 -4.03859901D3*TI 
     *         +3.81049965D0*TN(1) +4.06634913D-4*TN(2) 
     *         +5.21941168D-7*TN(3) -1.99565223D-10*TN(4) 
     *         +2.53447277D-14*TN(5) 
      ENDIF
C ho2cho
      IF (T .GT. 1378) THEN
      SMH(13) = -2.24939155D1 +3.80502496D4*TI 
     *         +9.87503878D0*TN(1) +2.32331854D-3*TN(2) 
     *         -2.78717537D-7*TN(3) +2.23853678D-11*TN(4) 
     *         -7.9797616D-16*TN(5) 
      ELSE
      SMH(13) = 1.75027796D1 +3.54828006D4*TI 
     *         +2.42464726D0*TN(1) +1.0985319D-2*TN(2) 
     *         -2.8117591D-6*TN(3) +5.21343495D-10*TN(4) 
     *         -4.55822921D-14*TN(5) 
      ENDIF
C o2cho
      IF (T .GT. 1368) THEN
      SMH(14) = -6.49547212D0 +1.87027618D4*TI 
     *         +7.24075139D0*TN(1) +2.31656475D-3*TN(2) 
     *         -2.72823325D-7*TN(3) +2.16422244D-11*TN(4) 
     *         -7.64823495D-16*TN(5) 
      ELSE
      SMH(14) = 1.17807483D1 +1.73599383D4*TI 
     *         +3.96059309D0*TN(1) +5.30011395D-3*TN(2) 
     *         -8.76188918D-7*TN(3) +8.47639383D-11*TN(4) 
     *         -1.43743801D-15*TN(5) 
      ENDIF
C ocho
      IF (T .GT. 1412) THEN
      SMH(15) = -1.00412613D1 +1.81679042D4*TI 
     *         +6.4981556D0*TN(1) +1.61159506D-3*TN(2) 
     *         -1.95165992D-7*TN(3) +1.57727432D-11*TN(4) 
     *         -5.64685D-16*TN(5) 
      ELSE
      SMH(15) = 1.80411779D1 +1.62042814D4*TI 
     *         +1.42991854D0*TN(1) +6.1401642D-3*TN(2) 
     *         -8.80768493D-7*TN(3) -3.10855459D-11*TN(4) 
     *         +2.53976786D-14*TN(5) 
      ENDIF
C ch2oh
      IF (T .GT. 1399) THEN
      SMH(16) = -3.49277963D0 +3.6147554D3*TI 
     *         +5.41875913D0*TN(1) +2.83092622D-3*TN(2) 
     *         -3.12451893D-7*TN(3) +2.37034948D-11*TN(4) 
     *         -8.1149899D-16*TN(5) 
      ELSE
      SMH(16) = 8.98878133D0 +2.8314019D3*TI 
     *         +3.05674228D0*TN(1) +5.9667817D-3*TN(2) 
     *         -1.45416884D-6*TN(3) +3.18984212D-10*TN(4) 
     *         -3.61443975D-14*TN(5) 
      ENDIF
C ch3o
      IF (T .GT. 1509) THEN
      SMH(17) = -1.57740193D0 +2.99208881D2*TI 
     *         +4.64787019D0*TN(1) +3.45415341D-3*TN(2) 
     *         -3.90674627D-7*TN(3) +3.01662142D-11*TN(4) 
     *         -1.0462677D-15*TN(5) 
      ELSE
      SMH(17) = 1.28377569D1 -9.45939708D2*TI 
     *         +2.23058023D0*TN(1) +4.26589293D-3*TN(2) 
     *         +1.70277707D-7*TN(3) -2.84205763D-10*TN(4) 
     *         +4.97345519D-14*TN(5) 
      ENDIF
C ch3o2h
      IF (T .GT. 1367) THEN
      SMH(18) = -2.17000591D1 +1.98512174D4*TI 
     *         +8.80409289D0*TN(1) +4.04713609D-3*TN(2) 
     *         -4.76405457D-7*TN(3) +3.77808128D-11*TN(4) 
     *         -1.33490354D-15*TN(5) 
      ELSE
      SMH(18) = 1.16092433D1 +1.74033753D4*TI 
     *         +2.83880024D0*TN(1) +9.30481245D-3*TN(2) 
     *         -1.41360902D-6*TN(3) +8.36562092D-11*TN(4) 
     *         +8.58062145D-15*TN(5) 
      ENDIF
C ch3o2
      IF (T .GT. 1365) THEN
      SMH(19) = -7.42552545D0 +1.83436055D3*TI 
     *         +6.34718801D0*TN(1) +3.96044679D-3*TN(2) 
     *         -4.61003188D-7*TN(3) +3.62800526D-11*TN(4) 
     *         -1.27492381D-15*TN(5) 
      ELSE
      SMH(19) = 7.817891D0 +4.55625796D2*TI 
     *         +3.8049759D0*TN(1) +4.9039233D-3*TN(2) 
     *         -6.51567707D-8*TN(3) -1.85893835D-10*TN(4) 
     *         +3.2165541D-14*TN(5) 
      ENDIF
C ch4
      IF (T .GT. 1462) THEN
      SMH(20) = -4.67561383D0 +1.13835704D4*TI 
     *         +4.09617653D0*TN(1) +3.72165422D-3*TN(2) 
     *         -4.397865D-7*TN(3) +3.49648003D-11*TN(4) 
     *         -1.23754025D-15*TN(5) 
      ELSE
      SMH(20) = 1.22776596D0 +1.01424099D4*TI 
     *         +3.7211302D0*TN(1) -1.25146645D-3*TN(2) 
     *         +3.17077557D-6*TN(3) -1.22392711D-9*TN(4) 
     *         +1.71895576D-13*TN(5) 
      ENDIF
C ch3
      IF (T .GT. 1389) THEN
      SMH(21) = 1.62436112D0 -1.61238027D4*TI 
     *         +3.51281376D0*TN(1) +2.55706306D-3*TN(2) 
     *         -2.7938675D-7*TN(3) +2.10412645D-11*TN(4) 
     *         -7.16514615D-16*TN(5) 
      ELSE
      SMH(21) = 2.52807406D0 -1.63164018D4*TI 
     *         +3.43858162D0*TN(1) +2.03876332D-3*TN(2) 
     *         +5.33051657D-8*TN(3) -7.89724492D-11*TN(4) 
     *         +1.10914083D-14*TN(5) 
      ENDIF
C ch2
      IF (T .GT. 1000) THEN
      SMH(22) = 2.156561D0 -4.534134D4*TI 
     *         +3.636408D0*TN(1) +9.665285D-4*TN(2) 
     *         -2.81169333D-8*TN(3) -8.415825D-12*TN(4) 
     *         +9.04128D-16*TN(5) 
      ELSE
      SMH(22) = 1.712578D0 -4.536791D4*TI 
     *         +3.762237D0*TN(1) +5.799095D-4*TN(2) 
     *         +4.14930833D-8*TN(3) +7.33403D-11*TN(4) 
     *         -3.6662175D-14*TN(5) 
      ENDIF
C ch2(s)
      IF (T .GT. 1000) THEN
      SMH(23) = 1.68657D0 -4.984975D4*TI 
     *         +3.552889D0*TN(1) +1.033394D-3*TN(2) 
     *         -3.19019333D-8*TN(3) -9.20560833D-12*TN(4) 
     *         +1.010675D-15*TN(5) 
      ELSE
      SMH(23) = 5.753207D-2 -4.989368D4*TI 
     *         +3.971265D0*TN(1) -8.495445D-5*TN(2) 
     *         +1.70894833D-7*TN(3) +2.07712583D-10*TN(4) 
     *         -9.90633D-14*TN(5) 
      ENDIF
C c2h5
      IF (T .GT. 1387) THEN
      SMH(24) = -8.49651771D0 -1.15065499D4*TI 
     *         +5.8878439D0*TN(1) +5.15383965D-3*TN(2) 
     *         -5.78073993D-7*TN(3) +4.43749381D-11*TN(4) 
     *         -1.53256325D-15*TN(5) 
      ELSE
      SMH(24) = 1.71789216D1 -1.34284028D4*TI 
     *         +1.32730217D0*TN(1) +8.83283765D-3*TN(2) 
     *         -1.0248776D-6*TN(3) -2.50952888D-11*TN(4) 
     *         +2.19308887D-14*TN(5) 
      ENDIF
C c2h4
      IF (T .GT. 1395) THEN
      SMH(25) = -7.47789234D0 -3.60389679D3*TI 
     *         +5.22176372D0*TN(1) +4.48068652D-3*TN(2) 
     *         -5.0811481D-7*TN(3) +3.92887937D-11*TN(4) 
     *         -1.36369796D-15*TN(5) 
      ELSE
      SMH(25) = 1.97084228D1 -5.46489338D3*TI 
     *         +2.33879687D-1*TN(1) +9.81673235D-3*TN(2) 
     *         -1.94722023D-6*TN(3) +3.03538711D-10*TN(4) 
     *         -2.38721358D-14*TN(5) 
      ENDIF
C c2h3
      IF (T .GT. 1395) THEN
      SMH(26) = -3.39792712D0 -3.37234748D4*TI 
     *         +5.07331248D0*TN(1) +3.29158139D-3*TN(2) 
     *         -3.72938207D-7*TN(3) +2.88169483D-11*TN(4) 
     *         -9.9970245D-16*TN(5) 
      ELSE
      SMH(26) = 1.71341661D1 -3.50734773D4*TI 
     *         +1.25329724D0*TN(1) +7.8129185D-3*TN(2) 
     *         -1.79673132D-6*TN(3) +3.48378862D-10*TN(4) 
     *         -3.50680181D-14*TN(5) 
      ENDIF
C ch3cho
      IF (T .GT. 1377) THEN
      SMH(27) = -1.27484852D1 +2.39807279D4*TI 
     *         +6.98518866D0*TN(1) +4.83948893D-3*TN(2) 
     *         -5.53069923D-7*TN(3) +4.30021584D-11*TN(4) 
     *         -1.49862952D-15*TN(5) 
      ELSE
      SMH(27) = 1.65023437D1 +2.1807885D4*TI 
     *         +1.77060035D0*TN(1) +9.22375805D-3*TN(2) 
     *         -1.20689694D-6*TN(3) +1.95303801D-11*TN(4) 
     *         +1.67771946D-14*TN(5) 
      ENDIF
C ch3co
      IF (T .GT. 1371) THEN
      SMH(28) = -8.83301965D0 +4.76690401D3*TI 
     *         +6.56682466D0*TN(1) +3.77654335D-3*TN(2) 
     *         -4.33277973D-7*TN(3) +3.37779162D-11*TN(4) 
     *         -1.17938082D-15*TN(5) 
      ELSE
      SMH(28) = 1.40340315D1 +3.02546532D3*TI 
     *         +2.5288415D0*TN(1) +6.85760865D-3*TN(2) 
     *         -7.14345793D-7*TN(3) -6.43070232D-11*TN(4) 
     *         +2.4191819D-14*TN(5) 
      ENDIF
C ch2cho
      IF (T .GT. 1388) THEN
      SMH(29) = -1.62602766D1 +2.61437239D3*TI 
     *         +7.5414579D0*TN(1) +3.41148562D-3*TN(2) 
     *         -3.94897678D-7*TN(3) +3.09694598D-11*TN(4) 
     *         -1.08580024D-15*TN(5) 
      ELSE
      SMH(29) = 1.63756465D1 +4.8205051D2*TI 
     *         +1.47616956D0*TN(1) +1.04487093D-2*TN(2) 
     *         -2.5020592D-6*TN(3) +4.69139485D-10*TN(4) 
     *         -4.38312235D-14*TN(5) 
      ENDIF
C c2h5oh
      IF (T .GT. 1395) THEN
      SMH(30) = -1.96502783D1 +3.25219715D4*TI 
     *         +8.31742137D0*TN(1) +6.4801647D-3*TN(2) 
     *         -7.32144162D-7*TN(3) +5.64779162D-11*TN(4) 
     *         -1.95724115D-15*TN(5) 
      ELSE
      SMH(30) = 2.44715919D1 +2.95670605D4*TI 
     *         +1.79106094D-1*TN(1) +1.54530041D-2*TN(2) 
     *         -3.22661162D-6*TN(3) +5.26525716D-10*TN(4) 
     *         -4.26583437D-14*TN(5) 
      ENDIF
C c2h5o
      IF (T .GT. 1393) THEN
      SMH(31) = -1.93190543D1 +6.22948597D3*TI 
     *         +8.23717244D0*TN(1) +5.54429395D-3*TN(2) 
     *         -6.31347287D-7*TN(3) +4.89677894D-11*TN(4) 
     *         -1.70356445D-15*TN(5) 
      ELSE
      SMH(31) = 2.37513898D1 +3.35717377D3*TI 
     *         +2.87429022D-1*TN(1) +1.43250459D-2*TN(2) 
     *         -3.06428347D-6*TN(3) +5.02580149D-10*TN(4) 
     *         -4.02281322D-14*TN(5) 
      ENDIF
C pc2h4oh
      IF (T .GT. 1396) THEN
      SMH(32) = -1.4467377D1 +7.44974622D3*TI 
     *         +7.88509591D0*TN(1) +5.41438975D-3*TN(2) 
     *         -6.11382805D-7*TN(3) +4.7144525D-11*TN(4) 
     *         -1.63336675D-15*TN(5) 
      ELSE
      SMH(32) = 2.28675716D1 +4.96032283D3*TI 
     *         +9.90023426D-1*TN(1) +1.30782136D-2*TN(2) 
     *         -2.76547167D-6*TN(3) +4.58770942D-10*TN(4) 
     *         -3.78760334D-14*TN(5) 
      ENDIF
C sc2h4oh
      IF (T .GT. 1397) THEN
      SMH(33) = -1.81354112D1 +1.10556188D4*TI 
     *         +8.34980768D0*TN(1) +5.1859034D-3*TN(2) 
     *         -5.82805802D-7*TN(3) +4.4793729D-11*TN(4) 
     *         -1.54830502D-15*TN(5) 
      ELSE
      SMH(33) = 1.80114775D1 +8.71494725D3*TI 
     *         +1.61103852D0*TN(1) +1.31755934D-2*TN(2) 
     *         -3.049025D-6*TN(3) +5.76035967D-10*TN(4) 
     *         -5.51920285D-14*TN(5) 
      ENDIF
C o2c2h4oh
      IF (T .GT. 1392) THEN
      SMH(34) = -2.33254953D1 +2.55911274D4*TI 
     *         +1.07432659D1*TN(1) +6.54788935D-3*TN(2) 
     *         -7.4228348D-7*TN(3) +5.73790615D-11*TN(4) 
     *         -1.99115056D-15*TN(5) 
      ELSE
      SMH(34) = 1.28482112D1 +2.30857785D4*TI 
     *         +4.11839445D0*TN(1) +1.36120316D-2*TN(2) 
     *         -2.68040717D-6*TN(3) +4.30861173D-10*TN(4) 
     *         -3.65805084D-14*TN(5) 
      ENDIF
C c2h3o1-2
      IF (T .GT. 1492) THEN
      SMH(35) = -1.2384257D1 -1.264422D4*TI 
     *         +6.88486471D0*TN(1) +3.4736025D-3*TN(2) 
     *         -3.72024497D-7*TN(3) +2.76825639D-11*TN(4) 
     *         -9.35622775D-16*TN(5) 
      ELSE
      SMH(35) = 3.22782741D1 -1.52459425D4*TI 
     *         -1.62965122D0*TN(1) +1.46727743D-2*TN(2) 
     *         -4.0622925D-6*TN(3) +8.37686042D-10*TN(4) 
     *         -8.0629518D-14*TN(5) 
      ENDIF
C n
      IF (T .GT. 1000) THEN
      SMH(36) = 4.448758D0 -5.611604D4*TI 
     *         +2.450268D0*TN(1) +5.33073D-5*TN(2) 
     *         -1.24422283D-8*TN(3) +1.56637667D-12*TN(4) 
     *         -5.12992D-17*TN(5) 
      ELSE
      SMH(36) = 4.167566D0 -5.60989D4*TI 
     *         +2.503071D0*TN(1) -1.090009D-5*TN(2) 
     *         +9.034215D-9*TN(3) -4.7063D-12*TN(4) 
     *         +1.049952D-15*TN(5) 
      ENDIF
C no
      IF (T .GT. 1000) THEN
      SMH(37) = 6.36900469D0 -9.89456954D3*TI 
     *         +3.26071234D0*TN(1) +5.95505675D-4*TN(2) 
     *         -7.1520441D-8*TN(3) +5.78734552D-12*TN(4) 
     *         -2.0164784D-16*TN(5) 
      ELSE
      SMH(37) = 2.28060952D0 -9.81823786D3*TI 
     *         +4.21859896D0*TN(1) -2.31994062D-3*TN(2) 
     *         +1.84071748D-6*TN(3) -7.78379589D-10*TN(4) 
     *         +1.40277437D-13*TN(5) 
      ENDIF
C no2
      IF (T .GT. 1000) THEN
      SMH(38) = -1.17416951D-1 -2.29397777D3*TI 
     *         +4.884754D0*TN(1) +1.08619775D-3*TN(2) 
     *         -1.38011515D-7*TN(3) +1.3122925D-11*TN(4) 
     *         -5.2554475D-16*TN(5) 
      ELSE
      SMH(38) = 6.3119919D0 -2.87409757D3*TI 
     *         +3.9440312D0*TN(1) -7.927145D-4*TN(2) 
     *         +2.776302D-6*TN(3) -1.7062855D-9*TN(4) 
     *         +3.9175282D-13*TN(5) 
      ENDIF
C n2o
      IF (T .GT. 1000) THEN
      SMH(39) = -1.65725D0 -8.165811D3*TI 
     *         +4.718977D0*TN(1) +1.436857D-3*TN(2) 
     *         -1.99582667D-7*TN(3) +1.87546D-11*TN(4) 
     *         -7.876685D-16*TN(5) 
      ELSE
      SMH(39) = 9.511222D0 -8.7651D3*TI 
     *         +2.543058D0*TN(1) +4.7460965D-3*TN(2) 
     *         -1.63212917D-6*TN(3) +5.21987083D-10*TN(4) 
     *         -9.50913D-14*TN(5) 
      ENDIF
C n2
      IF (T .GT. 1000) THEN
      SMH(40) = 5.980528D0 +9.227977D2*TI 
     *         +2.92664D0*TN(1) +7.439885D-4*TN(2) 
     *         -9.47460167D-8*TN(3) +8.4142D-12*TN(4) 
     *         -3.3766755D-16*TN(5) 
      ELSE
      SMH(40) = 3.950372D0 +1.0209D3*TI 
     *         +3.298677D0*TN(1) +7.0412D-4*TN(2) 
     *         -6.60537D-7*TN(3) +4.7012625D-10*TN(4) 
     *         -1.2224275D-13*TN(5) 
      ENDIF
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATX (T, C, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (SMALL = 1D-200)
      DIMENSION C(*), RF(*), RB(*), RKLOW(*)
C
      ALOGT = LOG(T)
      CTOT = 0.0
      DO K = 1, 28
         CTOT = CTOT + C(K)
      ENDDO
C
      RF(1) = RF(1)*C(1)*C(4)
      RB(1) = RB(1)*C(3)*C(5)
      RF(2) = RF(2)*C(2)*C(3)
      RB(2) = RB(2)*C(1)*C(5)
      RF(3) = RF(3)*C(2)*C(5)
      RB(3) = RB(3)*C(1)*C(6)
      RF(4) = RF(4)*C(3)*C(6)
      RB(4) = RB(4)*C(5)*C(5)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(6)+9D-1*C(9)+2.8D0*C(10)
      RF(5) = RF(5)*CTB*C(2)
      RB(5) = RB(5)*CTB*C(1)*C(1)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(6)+9D-1*C(9)+2.8D0*C(10)
     * +C(16)
      RF(6) = RF(6)*CTB*C(4)
      RB(6) = RB(6)*CTB*C(3)*C(3)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      RF(7) = RF(7)*CTB*C(5)
      RB(7) = RB(7)*CTB*C(1)*C(3)
      CTB = CTOT-2.7D-1*C(2)+1.1D1*C(6)+C(16)
      RF(8) = RF(8)*CTB*C(6)
      RB(8) = RB(8)*CTB*C(1)*C(5)
      CTB = CTOT+3D-1*C(2)+1.3D1*C(6)+9D-1*C(9)+2.8D0*C(10)
     * +C(16)
      PR = RKLOW(1) * CTB / RF(9)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 5D-1*EXP(-T/1D-30) + 5D-1*EXP(-T/1D30)
     *     + EXP(-1D30/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(9) = RF(9) * PCOR
      RB(9) = RB(9) * PCOR
      RF(9) = RF(9)*C(1)*C(4)
      RB(9) = RB(9)*C(7)
      RF(10) = RF(10)*C(1)*C(7)
      RB(10) = RB(10)*C(2)*C(4)
      RF(11) = RF(11)*C(1)*C(7)
      RB(11) = RB(11)*C(5)*C(5)
      RF(12) = RF(12)*C(3)*C(7)
      RB(12) = RB(12)*C(4)*C(5)
      RF(13) = RF(13)*C(5)*C(7)
      RB(13) = RB(13)*C(4)*C(6)
      RF(14) = RF(14)*C(4)*C(8)
      RB(14) = RB(14)*C(7)*C(7)
      RF(15) = RF(15)*C(4)*C(8)
      RB(15) = RB(15)*C(7)*C(7)
      CTB = CTOT+1.5D0*C(2)+1.1D1*C(6)+9D-1*C(9)+2.8D0*C(10)
     * +C(16)
      PR = RKLOW(2) * CTB / RF(16)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 5D-1*EXP(-T/1D-30) + 5D-1*EXP(-T/1D30)
     *     + EXP(-1D30/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(16) = RF(16) * PCOR
      RB(16) = RB(16) * PCOR
      RF(16) = RF(16)*C(8)
      RB(16) = RB(16)*C(5)*C(5)
      RF(17) = RF(17)*C(1)*C(8)
      RB(17) = RB(17)*C(5)*C(6)
      RF(18) = RF(18)*C(1)*C(8)
      RB(18) = RB(18)*C(2)*C(7)
      RF(19) = RF(19)*C(3)*C(8)
      RB(19) = RB(19)*C(5)*C(7)
      RF(20) = RF(20)*C(5)*C(8)
      RB(20) = RB(20)*C(6)*C(7)
      RF(21) = RF(21)*C(5)*C(8)
      RB(21) = RB(21)*C(6)*C(7)
      CTB = CTOT+C(2)+5D0*C(4)+5D0*C(6)+5D-1*C(9)
     * +2.5D0*C(10)+C(16)
      PR = RKLOW(3) * CTB / RF(22)
      PCOR = PR / (1.0 + PR)
      RF(22) = RF(22) * PCOR
      RB(22) = RB(22) * PCOR
      RF(22) = RF(22)*C(3)*C(9)
      RB(22) = RB(22)*C(10)
      RF(23) = RF(23)*C(4)*C(9)
      RB(23) = RB(23)*C(3)*C(10)
      RF(24) = RF(24)*C(5)*C(9)
      RB(24) = RB(24)*C(1)*C(10)
      RF(25) = RF(25)*C(7)*C(9)
      RB(25) = RB(25)*C(5)*C(10)
      CTB = CTOT+C(2)+1.1D1*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      RF(26) = RF(26)*CTB
      RB(26) = RB(26)*CTB*C(1)*C(9)
      RF(27) = RF(27)*C(4)
      RB(27) = RB(27)*C(7)*C(9)
      RF(28) = RF(28)*C(1)
      RB(28) = RB(28)*C(2)*C(9)
      RF(29) = RF(29)*C(3)
      RB(29) = RB(29)*C(5)*C(9)
      RF(30) = RF(30)*C(3)
      RB(30) = RB(30)*C(1)*C(10)
      RF(31) = RF(31)*C(5)
      RB(31) = RB(31)*C(6)*C(9)
      RF(32) = RF(32)*C(17)
      RB(32) = RB(32)*C(9)*C(16)
      RF(33) = RF(33)*C(7)
      RB(33) = RB(33)*C(4)*C(11)
      RF(34) = RF(34)*C(7)
      RB(34) = RB(34)*C(1)*C(5)*C(10)
      RF(35) = RF(35)*C(13)
      RB(35) = RB(35)*C(4)
      RF(36) = RF(36)*C(11)*C(13)
      RB(36) = RB(36)*C(12)
      RF(37) = RF(37)*C(12)
      RB(37) = RB(37)*C(5)
      RF(38) = RF(38)*C(1)*C(10)
      RF(39) = RF(39)*C(9)*C(11)
      RB(40) = RB(40)*C(2)*C(9)*C(9)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(4) * CTB / RF(41)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.176D-1*EXP(-T/2.71D2) + 7.824D-1*EXP(-T/2.755D3)
     *     + EXP(-6.57D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(41) = RF(41) * PCOR
      RB(41) = RB(41) * PCOR
      RF(41) = RF(41)*C(1)
      RB(41) = RB(41)*C(11)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(5) * CTB / RF(42)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 6.8D-2*EXP(-T/1.97D2) + 9.32D-1*EXP(-T/1.54D3)
     *     + EXP(-1.03D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(42) = RF(42) * PCOR
      RB(42) = RB(42) * PCOR
      RF(42) = RF(42)*C(2)*C(9)
      RB(42) = RB(42)*C(11)
      RF(43) = RF(43)*C(5)*C(11)
      RB(43) = RB(43)*C(6)
      RF(44) = RF(44)*C(1)*C(11)
      RB(44) = RB(44)*C(2)
      RF(45) = RF(45)*C(3)*C(11)
      RB(45) = RB(45)*C(5)
      RF(46) = RF(46)*C(11)*C(17)
      RB(46) = RB(46)*C(16)
      RF(47) = RF(47)*C(7)*C(11)
      RB(47) = RB(47)*C(8)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(6) * CTB / RF(48)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 10D-2*EXP(-T/2.5D3) + 9D-1*EXP(-T/1.3D3)
     *     + EXP(-10D98/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(48) = RF(48) * PCOR
      RB(48) = RB(48) * PCOR
      RB(48) = RB(48)*C(1)*C(11)
      RF(49) = RF(49)*C(4)
      RB(49) = RB(49)*C(7)*C(11)
      RF(50) = RF(50)*C(17)
      RB(50) = RB(50)*C(11)*C(16)
      RF(51) = RF(51)*C(1)
      RB(51) = RB(51)*C(2)*C(11)
      RF(52) = RF(52)*C(7)
      RB(52) = RB(52)*C(8)*C(11)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(7) * CTB / RF(53)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.813D-1*EXP(-T/1.03D2) + 7.187D-1*EXP(-T/1.291D3)
     *     + EXP(-4.16D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(53) = RF(53) * PCOR
      RB(53) = RB(53) * PCOR
      RF(53) = RF(53)*C(1)*C(11)
      RF(54) = RF(54)*C(4)
      RB(54) = RB(54)*C(7)*C(11)
      RF(55) = RF(55)*C(4)
      RB(55) = RB(55)*C(7)*C(11)
      RF(56) = RF(56)*C(1)
      RB(56) = RB(56)*C(2)*C(11)
      RF(57) = RF(57)*C(7)
      RB(57) = RB(57)*C(8)*C(11)
      RB(58) = RB(58)*C(11)*C(11)
      RF(59) = RF(59)*C(5)
      RB(59) = RB(59)*C(6)*C(11)
      RF(60) = RF(60)*C(3)
      RB(60) = RB(60)*C(5)*C(11)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(8) * CTB / RF(61)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.17D-1*EXP(-T/7.4D1) + 7.83D-1*EXP(-T/2.94D3)
     *     + EXP(-6.96D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(61) = RF(61) * PCOR
      RB(61) = RB(61) * PCOR
      RF(61) = RF(61)*C(1)*C(17)
      RB(61) = RB(61)*C(16)
      RF(62) = RF(62)*C(1)*C(16)
      RB(62) = RB(62)*C(2)*C(17)
      RF(63) = RF(63)*C(5)*C(16)
      RB(63) = RB(63)*C(6)*C(17)
      RF(64) = RF(64)*C(3)*C(16)
      RB(64) = RB(64)*C(5)*C(17)
      RF(65) = RF(65)*C(7)*C(16)
      RB(65) = RB(65)*C(8)*C(17)
      RF(66) = RF(66)*C(16)
      RB(66) = RB(66)*C(17)*C(17)
      RF(67) = RF(67)*C(5)*C(17)
      RB(67) = RB(67)*C(2)*C(11)
      RF(68) = RF(68)*C(5)*C(17)
      RB(68) = RB(68)*C(6)
      RF(69) = RF(69)*C(5)*C(17)
      RB(69) = RB(69)*C(1)
      RF(70) = RF(70)*C(5)*C(17)
      RB(70) = RB(70)*C(1)
      RF(71) = RF(71)*C(5)*C(17)
      RB(71) = RB(71)*C(6)
      RF(72) = RF(72)*C(7)*C(17)
      RB(72) = RB(72)*C(5)
      RF(73) = RF(73)*C(7)*C(17)
      RB(73) = RB(73)*C(4)*C(16)
      RF(74) = RF(74)*C(3)*C(17)
      RB(74) = RB(74)*C(1)*C(11)
      RF(75) = RF(75)*C(4)*C(17)
      RB(75) = RB(75)*C(3)
      RF(76) = RF(76)*C(4)*C(17)
      RB(76) = RB(76)*C(5)*C(11)
      CTB = CTOT
      PR = RKLOW(9) * CTB / RF(77)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4D-1*EXP(-T/1D3) + 6D-1*EXP(-T/7D1)
     *     + EXP(-1.7D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(77) = RF(77) * PCOR
      RB(77) = RB(77) * PCOR
      RF(77) = RF(77)*C(4)*C(17)
      RB(77) = RB(77)*C(15)
      RF(78) = RF(78)*C(11)*C(15)
      RB(78) = RB(78)*C(14)
      RF(79) = RF(79)*C(15)*C(16)
      RB(79) = RB(79)*C(14)*C(17)
      RF(80) = RF(80)*C(15)*C(17)
      RF(81) = RF(81)*C(7)*C(15)
      RB(81) = RB(81)*C(4)*C(14)
      RF(82) = RF(82)*C(15)*C(15)
      RB(82) = RB(82)*C(4)
      RF(83) = RF(83)*C(1)*C(15)
      RB(83) = RB(83)*C(5)
      RF(84) = RF(84)*C(3)*C(15)
      RB(84) = RB(84)*C(4)
      RF(85) = RF(85)*C(14)
      RB(85) = RB(85)*C(5)
      RF(87) = RF(87)*C(16)
      RB(87) = RB(87)*C(17)*C(17)
      RF(88) = RF(88)*C(4)
      RB(88) = RB(88)*C(1)*C(5)*C(9)
      RF(89) = RF(89)*C(2)
      RB(89) = RB(89)*C(1)*C(17)
      RF(90) = RF(90)*C(1)
      RB(90) = RB(90)*C(1)
      RF(91) = RF(91)*C(3)
      RB(91) = RB(91)*C(1)*C(1)*C(9)
      RF(92) = RF(92)*C(5)
      RB(92) = RB(92)*C(1)*C(11)
      RF(93) = RF(93)*C(10)
      RB(93) = RB(93)*C(9)*C(11)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(10) * CTB / RF(94)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.2D-1*EXP(-T/7.8D1) + 6.8D-1*EXP(-T/1.995D3)
     *     + EXP(-5.59D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(94) = RF(94) * PCOR
      RB(94) = RB(94) * PCOR
      RF(94) = RF(94)*C(1)
      RB(94) = RB(94)*C(17)
      RF(95) = RF(95)*C(4)
      RB(95) = RB(95)*C(3)*C(11)
      RF(96) = RF(96)*C(4)
      RB(96) = RB(96)*C(1)*C(1)*C(10)
      RF(97) = RF(97)*C(4)
      RB(97) = RB(97)*C(1)*C(5)*C(9)
      RF(98) = RF(98)*C(3)
      RB(98) = RB(98)*C(1)*C(1)*C(9)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(11) * CTB / RF(99)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.5D-2*EXP(-T/2.1D2) + 9.75D-1*EXP(-T/9.84D2)
     *     + EXP(-4.374D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(99) = RF(99) * PCOR
      RB(99) = RB(99) * PCOR
      RF(99) = RF(99)*C(1)*C(19)
      RB(99) = RB(99)*C(18)
      RF(100) = RF(100)*C(2)*C(15)
      RB(100) = RB(100)*C(1)*C(14)
      RF(101) = RF(101)*C(18)*C(20)
      RB(101) = RB(101)*C(19)*C(19)
      RF(102) = RF(102)*C(17)*C(18)
      RB(102) = RB(102)*C(16)*C(19)
      RF(103) = RF(103)*C(1)*C(18)
      RB(103) = RB(103)*C(17)*C(17)
      RF(104) = RF(104)*C(1)*C(18)
      RB(104) = RB(104)*C(2)*C(19)
      RF(105) = RF(105)*C(3)*C(18)
      RB(105) = RB(105)*C(1)*C(21)
      RF(106) = RF(106)*C(7)*C(18)
      RB(106) = RB(106)*C(5)
      RF(107) = RF(107)*C(15)*C(18)
      RF(108) = RF(108)*C(4)
      RB(108) = RB(108)*C(7)*C(21)
      RB(109) = RB(109)*C(11)*C(17)
      RB(110) = RB(110)*C(1)*C(21)
      RF(111) = RF(111)*C(4)*C(18)
      RB(111) = RB(111)*C(7)*C(19)
      RF(112) = RF(112)*C(4)*C(18)
      RB(112) = RB(112)*C(7)*C(19)
      RF(113) = RF(113)*C(4)*C(18)
      RB(113) = RB(113)*C(5)*C(21)
      RF(116) = RF(116)*C(21)
      RB(116) = RB(116)*C(17)
      RF(117) = RF(117)*C(1)*C(21)
      RB(117) = RB(117)*C(2)
      RF(118) = RF(118)*C(3)*C(21)
      RB(118) = RB(118)*C(5)
      RF(119) = RF(119)*C(5)*C(21)
      RB(119) = RB(119)*C(6)
      RF(120) = RF(120)*C(4)*C(21)
      RB(120) = RB(120)*C(7)
      RF(121) = RF(121)*C(17)*C(21)
      RB(121) = RB(121)*C(16)
      RF(122) = RF(122)*C(7)*C(21)
      RB(122) = RB(122)*C(8)
      RF(123) = RF(123)*C(15)*C(21)
      RB(123) = RB(123)*C(14)
      RF(124) = RF(124)*C(5)*C(21)
      RB(124) = RB(124)*C(6)
      CTB = CTOT
      PR = RKLOW(12) * CTB / RF(125)
      PCOR = PR / (1.0 + PR)
      RF(125) = RF(125) * PCOR
      RB(125) = RB(125) * PCOR
      RB(125) = RB(125)*C(9)*C(17)
      RF(126) = RF(126)*C(4)
      RB(126) = RB(126)*C(5)*C(9)*C(11)
      CTB = CTOT+C(2)+5D0*C(6)+5D-1*C(9)+C(10)
     * +C(16)
      PR = RKLOW(13) * CTB / RF(127)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.18D-1*EXP(-T/2.075D2) + 7.82D-1*EXP(-T/2.663D3)
     *     + EXP(-6.095D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(127) = RF(127) * PCOR
      RB(127) = RB(127) * PCOR
      RF(127) = RF(127)*C(1)*C(20)
      RB(127) = RB(127)*C(19)
      RF(128) = RF(128)*C(1)*C(19)
      RB(128) = RB(128)*C(2)*C(20)
      RF(129) = RF(129)*C(3)*C(19)
      RB(129) = RB(129)*C(17)
      RF(130) = RF(130)*C(3)*C(19)
      RB(130) = RB(130)*C(1)
      RF(131) = RF(131)*C(5)*C(19)
      RB(131) = RB(131)*C(6)*C(20)
      RF(132) = RF(132)*C(17)*C(19)
      RB(132) = RB(132)*C(16)*C(20)
      RF(133) = RF(133)*C(4)*C(19)
      RB(133) = RB(133)*C(7)*C(20)
      RF(134) = RF(134)*C(15)*C(19)
      RB(134) = RB(134)*C(14)*C(20)
      RF(135) = RF(135)*C(17)
      RB(135) = RB(135)*C(1)*C(19)
      RF(136) = RF(136)*C(4)*C(20)
      RB(136) = RB(136)*C(11)
      RF(137) = RF(137)*C(4)*C(20)
      RB(137) = RB(137)*C(3)
      CTB = CTOT+C(2)+4D0*C(6)+C(9)+2D0*C(10)
      PR = RKLOW(14) * CTB / RF(138)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 5D-1*EXP(-T/5.5D2) + 5D-1*EXP(-T/8.25D2)
     *     + EXP(-6.1D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(138) = RF(138) * PCOR
      RB(138) = RB(138) * PCOR
      RF(138) = RF(138)*C(22)
      RB(138) = RB(138)*C(17)
      CTB = CTOT+C(2)+4D0*C(6)+C(9)+2D0*C(10)
      PR = RKLOW(15) * CTB / RF(139)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 5D-1*EXP(-T/6.5D2) + 5D-1*EXP(-T/8D2)
     *     + EXP(-1D15/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(139) = RF(139) * PCOR
      RB(139) = RB(139) * PCOR
      RF(139) = RF(139)*C(22)
      RB(139) = RB(139)*C(5)*C(18)
      CTB = CTOT+4D0*C(6)
      PR = RKLOW(16) * CTB / RF(140)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 1.03D-1*EXP(-T/1D10) + 8.97D-1*EXP(-T/1D0)
     *     + EXP(-5D9/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(140) = RF(140) * PCOR
      RB(140) = RB(140) * PCOR
      RF(140) = RF(140)*C(22)
      RB(140) = RB(140)*C(6)*C(19)
      CTB = CTOT+4D0*C(6)
      PR = RKLOW(17) * CTB / RF(141)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 10D-2*EXP(-T/9D2) + 9D-1*EXP(-T/1.1D3)
     *     + EXP(-3.5D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(141) = RF(141) * PCOR
      RB(141) = RB(141) * PCOR
      RF(141) = RF(141)*C(22)
      RB(141) = RB(141)*C(2)*C(21)
      RF(142) = RF(142)*C(4)*C(22)
      RB(142) = RB(142)*C(7)
      RF(143) = RF(143)*C(4)*C(22)
      RB(143) = RB(143)*C(7)
      RF(144) = RF(144)*C(5)*C(22)
      RB(144) = RB(144)*C(6)
      RF(145) = RF(145)*C(5)*C(22)
      RB(145) = RB(145)*C(6)
      RF(146) = RF(146)*C(5)*C(22)
      RB(146) = RB(146)*C(6)
      RF(147) = RF(147)*C(1)*C(22)
      RB(147) = RB(147)*C(2)
      RF(148) = RF(148)*C(1)*C(22)
      RB(148) = RB(148)*C(2)
      RF(149) = RF(149)*C(1)*C(22)
      RB(149) = RB(149)*C(2)
      RF(150) = RF(150)*C(7)*C(22)
      RB(150) = RB(150)*C(8)
      RF(151) = RF(151)*C(7)*C(22)
      RB(151) = RB(151)*C(8)
      RF(152) = RF(152)*C(7)*C(22)
      RB(152) = RB(152)*C(8)
      RF(153) = RF(153)*C(15)*C(22)
      RB(153) = RB(153)*C(14)
      RF(154) = RF(154)*C(15)*C(22)
      RB(154) = RB(154)*C(14)
      RF(155) = RF(155)*C(15)*C(22)
      RB(155) = RB(155)*C(14)
      RF(156) = RF(156)*C(3)*C(22)
      RB(156) = RB(156)*C(5)
      RF(157) = RF(157)*C(3)*C(22)
      RB(157) = RB(157)*C(5)
      RF(158) = RF(158)*C(3)*C(22)
      RB(158) = RB(158)*C(5)
      RF(159) = RF(159)*C(17)*C(22)
      RB(159) = RB(159)*C(16)
      RF(160) = RF(160)*C(17)*C(22)
      RB(160) = RB(160)*C(16)
      RF(161) = RF(161)*C(17)*C(22)
      RB(161) = RB(161)*C(16)
      RB(162) = RB(162)*C(5)*C(19)
      CTB = CTOT
      RF(163) = RF(163)*CTB
      RB(163) = RB(163)*CTB*C(1)*C(21)
      RF(164) = RF(164)*C(23)
      RB(164) = RB(164)*C(4)
      RF(165) = RF(165)*C(23)
      RB(165) = RB(165)*C(5)*C(11)*C(11)
      RF(166) = RF(166)*C(4)
      RB(166) = RB(166)*C(7)*C(21)
      RF(167) = RF(167)*C(5)*C(24)
      RB(167) = RB(167)*C(1)*C(25)
      RF(168) = RF(168)*C(4)*C(24)
      RB(168) = RB(168)*C(3)*C(25)
      RF(169) = RF(169)*C(24)*C(25)
      RB(169) = RB(169)*C(3)*C(28)
      CTB = CTOT
      PR = RKLOW(18) * CTB / RF(170)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 1.2D-1*EXP(-T/1D3) + 8.8D-1*EXP(-T/1D4)
     *     + EXP(-1D30/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(170) = RF(170) * PCOR
      RB(170) = RB(170) * PCOR
      RF(170) = RF(170)*C(3)*C(25)
      RB(170) = RB(170)*C(26)
      RF(171) = RF(171)*C(7)*C(25)
      RB(171) = RB(171)*C(5)*C(26)
      RF(172) = RF(172)*C(1)*C(26)
      RB(172) = RB(172)*C(5)*C(25)
      RF(173) = RF(173)*C(3)*C(26)
      RB(173) = RB(173)*C(4)*C(25)
      RF(174) = RF(174)*C(26)*C(26)
      RB(174) = RB(174)*C(4)*C(25)*C(25)
      CTB = CTOT+7D-1*C(28)+4D-1*C(4)+2D0*C(10)+1.1D1*C(6)
      PR = RKLOW(19) * CTB / RF(175)
      PCOR = PR / (1.0 + PR)
      RF(175) = RF(175) * PCOR
      RB(175) = RB(175) * PCOR
      RF(175) = RF(175)*C(27)
      RB(175) = RB(175)*C(3)*C(28)
      RF(176) = RF(176)*C(1)*C(27)
      RB(176) = RB(176)*C(5)*C(28)
      RF(177) = RF(177)*C(3)*C(27)
      RB(177) = RB(177)*C(25)*C(25)
      RF(178) = RF(178)*C(3)*C(27)
      RB(178) = RB(178)*C(4)*C(28)
      RF(179) = RF(179)*C(5)*C(27)
      RB(179) = RB(179)*C(7)*C(28)
      RF(180) = RF(180)*C(25)*C(27)
      RB(180) = RB(180)*C(26)*C(28)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDOT(RF, RB, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), WDOT(*)
C
      DO K = 1, 28
         WDOT(K) = 0D0
      ENDDO
C
      ROP = RF(1)-RB(1)
      WDOT(1) = WDOT(1) -ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      ROP = RF(2)-RB(2)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      ROP = RF(3)-RB(3)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(4)-RB(4)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP +ROP
      WDOT(6) = WDOT(6) -ROP
      ROP = RF(5)-RB(5)
      WDOT(1) = WDOT(1) +ROP +ROP
      WDOT(2) = WDOT(2) -ROP
      ROP = RF(6)-RB(6)
      WDOT(3) = WDOT(3) +ROP +ROP
      WDOT(4) = WDOT(4) -ROP
      ROP = RF(7)-RB(7)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(5) = WDOT(5) -ROP
      ROP = RF(8)-RB(8)
      WDOT(1) = WDOT(1) +ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(6) = WDOT(6) -ROP
      ROP = RF(9)-RB(9)
      WDOT(1) = WDOT(1) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(10)-RB(10)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(11)-RB(11)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(12)-RB(12)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(13)-RB(13)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(14)-RB(14)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(15)-RB(15)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(16)-RB(16)
      WDOT(5) = WDOT(5) +ROP +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(17)-RB(17)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(18)-RB(18)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(19)-RB(19)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(20)-RB(20)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(21)-RB(21)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(22)-RB(22)
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(23)-RB(23)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(24)-RB(24)
      WDOT(1) = WDOT(1) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(25)-RB(25)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(26)-RB(26)
      WDOT(1) = WDOT(1) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(27)-RB(27)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(28)-RB(28)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(29)-RB(29)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(30)-RB(30)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(31)-RB(31)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(32)-RB(32)
      WDOT(9) = WDOT(9) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(33)-RB(33)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(34)-RB(34)
      WDOT(1) = WDOT(1) +ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(35)-RB(35)
      WDOT(4) = WDOT(4) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(36)-RB(36)
      WDOT(11) = WDOT(11) -ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(37)-RB(37)
      WDOT(5) = WDOT(5) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(38)-RB(38)
      WDOT(1) = WDOT(1) -ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(39)-RB(39)
      WDOT(9) = WDOT(9) -ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(40)-RB(40)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) +ROP +ROP
      ROP = RF(41)-RB(41)
      WDOT(1) = WDOT(1) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(42)-RB(42)
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(43)-RB(43)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(44)-RB(44)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(45)-RB(45)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(46)-RB(46)
      WDOT(11) = WDOT(11) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(47)-RB(47)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(48)-RB(48)
      WDOT(1) = WDOT(1) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(49)-RB(49)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(50)-RB(50)
      WDOT(11) = WDOT(11) +ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(51)-RB(51)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(52)-RB(52)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(53)-RB(53)
      WDOT(1) = WDOT(1) -ROP
      WDOT(11) = WDOT(11) -ROP
      ROP = RF(54)-RB(54)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(55)-RB(55)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(56)-RB(56)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(57)-RB(57)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(58)-RB(58)
      WDOT(11) = WDOT(11) +ROP +ROP
      ROP = RF(59)-RB(59)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(60)-RB(60)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(61)-RB(61)
      WDOT(1) = WDOT(1) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(62)-RB(62)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(63)-RB(63)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(64)-RB(64)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(65)-RB(65)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(66)-RB(66)
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP +ROP
      ROP = RF(67)-RB(67)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(68)-RB(68)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(69)-RB(69)
      WDOT(1) = WDOT(1) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(70)-RB(70)
      WDOT(1) = WDOT(1) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(71)-RB(71)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(72)-RB(72)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(73)-RB(73)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(74)-RB(74)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(75)-RB(75)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(76)-RB(76)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(77)-RB(77)
      WDOT(4) = WDOT(4) -ROP
      WDOT(15) = WDOT(15) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(78)-RB(78)
      WDOT(11) = WDOT(11) -ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(79)-RB(79)
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(80)-RB(80)
      WDOT(15) = WDOT(15) -ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(81)-RB(81)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(82)-RB(82)
      WDOT(4) = WDOT(4) +ROP
      WDOT(15) = WDOT(15) -ROP -ROP
      ROP = RF(83)-RB(83)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(84)-RB(84)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(85)-RB(85)
      WDOT(5) = WDOT(5) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(86)-RB(86)
      ROP = RF(87)-RB(87)
      WDOT(16) = WDOT(16) -ROP
      WDOT(17) = WDOT(17) +ROP +ROP
      ROP = RF(88)-RB(88)
      WDOT(1) = WDOT(1) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(89)-RB(89)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(90)-RB(90)
      ROP = RF(91)-RB(91)
      WDOT(1) = WDOT(1) +ROP +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(92)-RB(92)
      WDOT(1) = WDOT(1) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(93)-RB(93)
      WDOT(9) = WDOT(9) +ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(94)-RB(94)
      WDOT(1) = WDOT(1) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(95)-RB(95)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(96)-RB(96)
      WDOT(1) = WDOT(1) +ROP +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(97)-RB(97)
      WDOT(1) = WDOT(1) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(98)-RB(98)
      WDOT(1) = WDOT(1) +ROP +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(99)-RB(99)
      WDOT(1) = WDOT(1) -ROP
      WDOT(18) = WDOT(18) +ROP
      WDOT(19) = WDOT(19) -ROP
      ROP = RF(100)-RB(100)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(101)-RB(101)
      WDOT(18) = WDOT(18) -ROP
      WDOT(19) = WDOT(19) +ROP +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(102)-RB(102)
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      WDOT(18) = WDOT(18) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(103)-RB(103)
      WDOT(1) = WDOT(1) -ROP
      WDOT(17) = WDOT(17) +ROP +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(104)-RB(104)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(18) = WDOT(18) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(105)-RB(105)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(18) = WDOT(18) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(106)-RB(106)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(107)-RB(107)
      WDOT(15) = WDOT(15) -ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(108)-RB(108)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(109)-RB(109)
      WDOT(11) = WDOT(11) +ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(110)-RB(110)
      WDOT(1) = WDOT(1) +ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(111)-RB(111)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(18) = WDOT(18) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(112)-RB(112)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(18) = WDOT(18) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(113)-RB(113)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(18) = WDOT(18) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(114)-RB(114)
      ROP = RF(115)-RB(115)
      ROP = RF(116)-RB(116)
      WDOT(17) = WDOT(17) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(117)-RB(117)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(118)-RB(118)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(119)-RB(119)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(120)-RB(120)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(121)-RB(121)
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(122)-RB(122)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(123)-RB(123)
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(124)-RB(124)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(125)-RB(125)
      WDOT(9) = WDOT(9) +ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(126)-RB(126)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(127)-RB(127)
      WDOT(1) = WDOT(1) -ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(128)-RB(128)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(129)-RB(129)
      WDOT(3) = WDOT(3) -ROP
      WDOT(17) = WDOT(17) +ROP
      WDOT(19) = WDOT(19) -ROP
      ROP = RF(130)-RB(130)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(19) = WDOT(19) -ROP
      ROP = RF(131)-RB(131)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(132)-RB(132)
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(133)-RB(133)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(134)-RB(134)
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(19) = WDOT(19) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(135)-RB(135)
      WDOT(1) = WDOT(1) +ROP
      WDOT(17) = WDOT(17) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(136)-RB(136)
      WDOT(4) = WDOT(4) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(137)-RB(137)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(138)-RB(138)
      WDOT(17) = WDOT(17) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(139)-RB(139)
      WDOT(5) = WDOT(5) +ROP
      WDOT(18) = WDOT(18) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(140)-RB(140)
      WDOT(6) = WDOT(6) +ROP
      WDOT(19) = WDOT(19) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(141)-RB(141)
      WDOT(2) = WDOT(2) +ROP
      WDOT(21) = WDOT(21) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(142)-RB(142)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(143)-RB(143)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(144)-RB(144)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(145)-RB(145)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(146)-RB(146)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(147)-RB(147)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(148)-RB(148)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(149)-RB(149)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(150)-RB(150)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(151)-RB(151)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(152)-RB(152)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(153)-RB(153)
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(154)-RB(154)
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(155)-RB(155)
      WDOT(14) = WDOT(14) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(156)-RB(156)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(157)-RB(157)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(158)-RB(158)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(159)-RB(159)
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(160)-RB(160)
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(161)-RB(161)
      WDOT(16) = WDOT(16) +ROP
      WDOT(17) = WDOT(17) -ROP
      WDOT(22) = WDOT(22) -ROP
      ROP = RF(162)-RB(162)
      WDOT(5) = WDOT(5) +ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(163)-RB(163)
      WDOT(1) = WDOT(1) +ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(164)-RB(164)
      WDOT(4) = WDOT(4) +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(165)-RB(165)
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP +ROP
      WDOT(23) = WDOT(23) -ROP
      ROP = RF(166)-RB(166)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(167)-RB(167)
      WDOT(1) = WDOT(1) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(24) = WDOT(24) -ROP
      WDOT(25) = WDOT(25) +ROP
      ROP = RF(168)-RB(168)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(24) = WDOT(24) -ROP
      WDOT(25) = WDOT(25) +ROP
      ROP = RF(169)-RB(169)
      WDOT(3) = WDOT(3) +ROP
      WDOT(24) = WDOT(24) -ROP
      WDOT(25) = WDOT(25) -ROP
      WDOT(28) = WDOT(28) +ROP
      ROP = RF(170)-RB(170)
      WDOT(3) = WDOT(3) -ROP
      WDOT(25) = WDOT(25) -ROP
      WDOT(26) = WDOT(26) +ROP
      ROP = RF(171)-RB(171)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(25) = WDOT(25) -ROP
      WDOT(26) = WDOT(26) +ROP
      ROP = RF(172)-RB(172)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(25) = WDOT(25) +ROP
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(173)-RB(173)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(25) = WDOT(25) +ROP
      WDOT(26) = WDOT(26) -ROP
      ROP = RF(174)-RB(174)
      WDOT(4) = WDOT(4) +ROP
      WDOT(25) = WDOT(25) +ROP +ROP
      WDOT(26) = WDOT(26) -ROP -ROP
      ROP = RF(175)-RB(175)
      WDOT(3) = WDOT(3) +ROP
      WDOT(27) = WDOT(27) -ROP
      WDOT(28) = WDOT(28) +ROP
      ROP = RF(176)-RB(176)
      WDOT(1) = WDOT(1) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(27) = WDOT(27) -ROP
      WDOT(28) = WDOT(28) +ROP
      ROP = RF(177)-RB(177)
      WDOT(3) = WDOT(3) -ROP
      WDOT(25) = WDOT(25) +ROP +ROP
      WDOT(27) = WDOT(27) -ROP
      ROP = RF(178)-RB(178)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(27) = WDOT(27) -ROP
      WDOT(28) = WDOT(28) +ROP
      ROP = RF(179)-RB(179)
      WDOT(5) = WDOT(5) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(27) = WDOT(27) -ROP
      WDOT(28) = WDOT(28) +ROP
      ROP = RF(180)-RB(180)
      WDOT(25) = WDOT(25) -ROP
      WDOT(26) = WDOT(26) +ROP
      WDOT(27) = WDOT(27) -ROP
      WDOT(28) = WDOT(28) +ROP
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE QSSA(RF, RB, XQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), XQ(*)
C
      RB(39) = 0.D0
      RF(40) = 0.D0
      RF(58) = 0.D0
      RB(80) = 0.D0
      RB(82) = 0.D0
      RB(107) = 0.D0
C
C     hco
      DEN = +RF( 26) +RF( 27) +RF( 28) +RF( 29) +RF( 30) 
     *  +RF( 31) +RF( 32) +RF( 33) +RF( 34) +RF( 41) +RB( 35) 
     *  +RB( 36) +RB( 43) +RB( 44) +RB( 45) +RB( 46) +RB( 47) 
     *  +RB( 78) +RB(116) +RB(129) +RB(136) 
      A1_0 = ( +RB( 26) +RB( 27) +RB( 28) +RB( 29) +RB( 30) 
     *  +RB( 31) +RB( 32) +RB( 33) +RB( 34) +RF( 35) +RF( 36) 
     *  +RF( 39) +RF( 39) +RB( 40) +RB( 40) +RB( 41) +RF( 43) 
     *  +RF( 44) +RF( 45) +RF( 46) +RF( 47) +RB( 58) +RF( 78) 
     *  +RF(116) +RF(129) +RF(136) )/DEN
C     ocho
      DEN = +RB( 37) +RB( 38) 
      A2_0 = ( +RF( 37) +RF( 38) )/DEN
C     ch2oh
      DEN = +RF( 54) +RF( 55) +RF( 56) +RF( 57) +RF( 59) 
     *  +RF( 60) +RB( 53) +RB( 70) +RB(138) 
      A3_0 = ( +RF( 53) +RB( 54) +RB( 55) +RB( 56) +RB( 57) 
     *  +RB( 58) +RB( 59) +RB( 60) +RF( 70) +RF(138) )/DEN
C     ch3o
      DEN = +RF( 48) +RF( 49) +RF( 50) +RF( 51) +RF( 52) 
     *  +RB( 69) +RB( 72) +RB( 75) +RB( 83) +RB( 84) +RB( 85) 
      A4_0 = ( +RB( 48) +RB( 49) +RB( 50) +RB( 51) +RB( 52) 
     *  +RF( 69) +RF( 72) +RF( 75) +RF( 80) +RF( 80) +RF( 82) 
     *  +RF( 82) +RF( 83) +RF( 84) +RF( 85) +RF(107) )/DEN
C     ch2
      DEN = +RF( 66) +RF( 94) +RF( 95) +RF( 96) +RF( 97) 
     *  +RF( 98) +RB( 71) +RB( 86) +RB( 90) 
      A5_0 = ( +RB( 66) +RF( 71) +RB( 94) +RB( 95) +RB( 96) 
     *  +RB( 97) +RB( 98) )/DEN
      A5_6 = ( +RF( 86) +RF( 90) )/DEN
C     ch2(s)
      DEN = +RF( 86) +RF( 87) +RF( 88) +RF( 89) +RF( 90) 
     *  +RF( 91) +RF( 92) +RF( 93) +RF(135) +RB( 68) 
      A6_0 = ( +RF( 68) +RB( 87) +RB( 88) +RB( 89) +RB( 91) 
     *  +RB( 92) +RB( 93) +RB(135) )/DEN
      A6_5 = ( +RB( 86) +RB( 90) )/DEN
C     ch3co
      DEN = +RF(125) +RB(114) +RB(117) +RB(118) +RB(119) 
     *  +RB(120) +RB(121) +RB(122) +RB(123) 
      A7_0 = ( +RF(117) +RF(118) +RF(119) +RF(120) +RF(121) 
     *  +RF(122) +RF(123) +RB(125) )/DEN
      A7_12 = ( +RF(114) )/DEN
C     ch2cho
      DEN = +RF(126) +RB(115) +RB(124) +RB(130) +RB(137) 
      A8_0 = ( +RF(124) +RB(126) +RF(130) +RF(137) )/DEN
      A8_12 = ( +RF(115) )/DEN
C     c2h5o
      DEN = +RF(108) +RF(109) +RF(110) +RB(106) +RB(146) 
     *  +RB(149) +RB(152) +RB(155) +RB(158) +RB(161) 
      A9_0 = ( +RF(106) +RF(107) +RB(108) +RB(109) +RB(110) 
     *  +RF(146) +RF(149) +RF(152) +RF(155) +RF(158) +RF(161) )/DEN
C     pc2h4oh
      DEN = +RF(162) +RB(142) +RB(144) +RB(147) +RB(150) 
     *  +RB(153) +RB(156) +RB(159) +RB(164) 
      A10_0 = ( +RF(142) +RF(144) +RF(147) +RF(150) +RF(153) 
     *  +RF(156) +RF(159) +RB(162) +RF(164) )/DEN
C     sc2h4oh
      DEN = +RF(163) +RF(166) +RB(143) +RB(145) +RB(148) 
     *  +RB(151) +RB(154) +RB(157) +RB(160) 
      A11_0 = ( +RF(143) +RF(145) +RF(148) +RF(151) +RF(154) 
     *  +RF(157) +RF(160) +RB(163) +RB(166) )/DEN
C     c2h3o1-2
      DEN = +RF(114) +RF(115) 
C      A12_0 = ( )/DEN
      A12_0 = 0.D0
      A12_7 = ( +RB(114) )/DEN
      A12_8 = ( +RB(115) )/DEN
C
      XQ(1) = A1_0
      XQ(2) = A2_0
      XQ(3) = A3_0
      XQ(4) = A4_0
      A5_0 = A5_0 + A5_6*A6_0
      DEN = 1 -A5_6*A6_5
      A5_0 = A5_0/DEN
      XQ(5) = A5_0
      XQ(6) = A6_0 +A6_5*XQ(5)
      A12_0 = A12_0 + A12_7*A7_0
      DEN = 1 -A12_7*A7_12
      A12_0 = A12_0/DEN
      A12_8 = A12_8/DEN
      A12_0 = A12_0 + A12_8*A8_0
      DEN = 1 -A12_8*A8_12
      A12_0 = A12_0/DEN
      XQ(12) = A12_0
      XQ(8) = A8_0 +A8_12*XQ(12)
      XQ(7) = A7_0 +A7_12*XQ(12)
      XQ(9) = A9_0
      XQ(10) = A10_0
      XQ(11) = A11_0
C
      RF( 26) = RF( 26)*XQ( 1)
      RF( 27) = RF( 27)*XQ( 1)
      RF( 28) = RF( 28)*XQ( 1)
      RF( 29) = RF( 29)*XQ( 1)
      RF( 30) = RF( 30)*XQ( 1)
      RF( 31) = RF( 31)*XQ( 1)
      RF( 32) = RF( 32)*XQ( 1)
      RF( 33) = RF( 33)*XQ( 1)
      RF( 34) = RF( 34)*XQ( 1)
      RB( 35) = RB( 35)*XQ( 1)
      RB( 36) = RB( 36)*XQ( 1)
      RB( 37) = RB( 37)*XQ( 2)
      RB( 38) = RB( 38)*XQ( 2)
      RF( 41) = RF( 41)*XQ( 1)
      RB( 43) = RB( 43)*XQ( 1)
      RB( 44) = RB( 44)*XQ( 1)
      RB( 45) = RB( 45)*XQ( 1)
      RB( 46) = RB( 46)*XQ( 1)
      RB( 47) = RB( 47)*XQ( 1)
      RF( 48) = RF( 48)*XQ( 4)
      RF( 49) = RF( 49)*XQ( 4)
      RF( 50) = RF( 50)*XQ( 4)
      RF( 51) = RF( 51)*XQ( 4)
      RF( 52) = RF( 52)*XQ( 4)
      RB( 53) = RB( 53)*XQ( 3)
      RF( 54) = RF( 54)*XQ( 3)
      RF( 55) = RF( 55)*XQ( 3)
      RF( 56) = RF( 56)*XQ( 3)
      RF( 57) = RF( 57)*XQ( 3)
      RF( 59) = RF( 59)*XQ( 3)
      RF( 60) = RF( 60)*XQ( 3)
      RF( 66) = RF( 66)*XQ( 5)
      RB( 68) = RB( 68)*XQ( 6)
      RB( 69) = RB( 69)*XQ( 4)
      RB( 70) = RB( 70)*XQ( 3)
      RB( 71) = RB( 71)*XQ( 5)
      RB( 72) = RB( 72)*XQ( 4)
      RB( 75) = RB( 75)*XQ( 4)
      RB( 78) = RB( 78)*XQ( 1)
      RB( 83) = RB( 83)*XQ( 4)
      RB( 84) = RB( 84)*XQ( 4)
      RB( 85) = RB( 85)*XQ( 4)
      RF( 86) = RF( 86)*XQ( 6)
      RB( 86) = RB( 86)*XQ( 5)
      RF( 87) = RF( 87)*XQ( 6)
      RF( 88) = RF( 88)*XQ( 6)
      RF( 89) = RF( 89)*XQ( 6)
      RF( 90) = RF( 90)*XQ( 6)
      RB( 90) = RB( 90)*XQ( 5)
      RF( 91) = RF( 91)*XQ( 6)
      RF( 92) = RF( 92)*XQ( 6)
      RF( 93) = RF( 93)*XQ( 6)
      RF( 94) = RF( 94)*XQ( 5)
      RF( 95) = RF( 95)*XQ( 5)
      RF( 96) = RF( 96)*XQ( 5)
      RF( 97) = RF( 97)*XQ( 5)
      RF( 98) = RF( 98)*XQ( 5)
      RB(106) = RB(106)*XQ( 9)
      RF(108) = RF(108)*XQ( 9)
      RF(109) = RF(109)*XQ( 9)
      RF(110) = RF(110)*XQ( 9)
      RF(114) = RF(114)*XQ(12)
      RB(114) = RB(114)*XQ( 7)
      RF(115) = RF(115)*XQ(12)
      RB(115) = RB(115)*XQ( 8)
      RB(116) = RB(116)*XQ( 1)
      RB(117) = RB(117)*XQ( 7)
      RB(118) = RB(118)*XQ( 7)
      RB(119) = RB(119)*XQ( 7)
      RB(120) = RB(120)*XQ( 7)
      RB(121) = RB(121)*XQ( 7)
      RB(122) = RB(122)*XQ( 7)
      RB(123) = RB(123)*XQ( 7)
      RB(124) = RB(124)*XQ( 8)
      RF(125) = RF(125)*XQ( 7)
      RF(126) = RF(126)*XQ( 8)
      RB(129) = RB(129)*XQ( 1)
      RB(130) = RB(130)*XQ( 8)
      RF(135) = RF(135)*XQ( 6)
      RB(136) = RB(136)*XQ( 1)
      RB(137) = RB(137)*XQ( 8)
      RB(138) = RB(138)*XQ( 3)
      RB(142) = RB(142)*XQ(10)
      RB(143) = RB(143)*XQ(11)
      RB(144) = RB(144)*XQ(10)
      RB(145) = RB(145)*XQ(11)
      RB(146) = RB(146)*XQ( 9)
      RB(147) = RB(147)*XQ(10)
      RB(148) = RB(148)*XQ(11)
      RB(149) = RB(149)*XQ( 9)
      RB(150) = RB(150)*XQ(10)
      RB(151) = RB(151)*XQ(11)
      RB(152) = RB(152)*XQ( 9)
      RB(153) = RB(153)*XQ(10)
      RB(154) = RB(154)*XQ(11)
      RB(155) = RB(155)*XQ( 9)
      RB(156) = RB(156)*XQ(10)
      RB(157) = RB(157)*XQ(11)
      RB(158) = RB(158)*XQ( 9)
      RB(159) = RB(159)*XQ(10)
      RB(160) = RB(160)*XQ(11)
      RB(161) = RB(161)*XQ( 9)
      RF(162) = RF(162)*XQ(10)
      RF(163) = RF(163)*XQ(11)
      RB(164) = RB(164)*XQ(10)
      RF(166) = RF(166)*XQ(11)
C
      END
