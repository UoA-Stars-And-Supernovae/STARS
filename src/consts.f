      SUBROUTINE CONSTS
      IMPLICIT REAL*8 (A-H,L,M,O-Z)
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CL, CD, CG, CR, CT, CEVB, 
     &     CEN, CPL, CMEVMU, CSECYR, LSUN, MSUN, RSUN, TSUNYR,
     &     STEFBOLTZ
      COMMON /ATDATA/ DH2(4), CHI(26,9), OMG(27), AM(10), BN(10), IZ(10)
      COMMON /NCDATA/ QRT(20), QNT(20), CZA(92), CZB(92), CZC(92),
     &     CZD(92), VZ(10)
      COMMON /MASLOS/ AIJ(6,5),baseN
      DIMENSION Z1(92),Z2(92)
* original values for constants
      DATA CPI, PI4, CLN10/3.1416D0, 12.5664D0, 2.3026D0/
      DATA CA, CB, CL, CD, CG, CR, CT/7.5647E-15, 1.4406E24, 2.9979D10,
     &     2.9218E6, 6.672D-8, 8.3143E7, 1.6863E-10/
      DATA CEVB, CSECYR, CMEVMU, LSUN, MSUN, RSUN, TSUNYR /1.1605D4,
     &     3.1557D7, 9.6488D17, 3.8515D0, 1.9891D0, 0.69598D0, 4.57D9/
      DATA STEFBOLTZ /5.67037D-8/
* mathematical constants
      CPI = 4*ATAN(1.D0)
      PI4 = 4*CPI
      CLN10 = LOG(10.0D0)
* physical constants in cgs units, copied from JCD's code
      AMU = 1.6605402D-24
      AME = 9.1093897D-28
      CL = 2.99792458D10
      PLANCK = 6.6260755D-27
      BOLTZM = 1.380658D-16
      ECHAR = 4.8032068D-10
      EVOLT = 1.60217733D-12
      CSECYR = 3.155692597D7
      CG = 6.672D-8
* derived physical constants
      CA = BOLTZM/(CL*PLANCK)
      CA = 8*CPI**5*BOLTZM*CA**3/15
      CMEVMU = 1.0D6*EVOLT/AMU
      CEVB = EVOLT/BOLTZM
      CR = BOLTZM/AMU
      CHE = AME*CL**2
      CT = BOLTZM/CHE
      LAMC = PLANCK/(AME*CL)
      LAMC3 = LAMC**3
      CRHO = 8*CPI/LAMC3
      CB = CRHO*CHE
      CD = CRHO*AMU
      CEN = PLANCK/SQRT(2*CPI*AMU*BOLTZM)
      CEN = CEN**3/AMU
      CPL = SQRT(PI4/(AMU*BOLTZM**3))*ECHAR**3
      STEFBOLTZ = 2d0 * CPI**5d0 * BOLTZM**4 / (15d0 * PLANCK**3d0 * CL**2)
* solar mass, luminosity, radius, age
      MSUN = 1.9891D0
      LSUN = 3.844D0
      RSUN = 0.69598D0
      TSUNYR = 4.57D9
* nuclear reaction and neutrino Q values, in MeV
C      DATA QRT/6.936, 12.86, 1.586, 18.21, 18.21, 11.715, 15.017, 4.553,
      DATA QRT/6.936, 12.86, 19.796, 18.21, 18.21, 11.715, 15.017, 4.553,
     &         7.275, 7.162, 4.415, 4.734, 9.312, 4.261, 6.771, -0.391,
     &         -4.734, -9.312, 13.933, 22.179/
C      DATA QNT/0.265, 0.0, 0.0, 0.814, 6.710, 0.707, 0.997, 0.999,
      DATA QNT/0.265, 0.0, 0.814, 0.814, 6.710, 0.707, 0.997, 0.999,
     &         11*0.0, 0.997/
* constants for electron screening
      DATA Z1 /1,2,2,4,4,6,7,8,4,6,7,8,10,6,8,8,10,12,0,0,1,1,1,0,1,1,1,
     :     2,2,2,2,1,1,2,1,1,2,2,1,1,2,1,2,2,1,2,2,0,0,1,1,1,1,1,4,1,1,1,
     :     4,4,4,1,1,1,1,4,4,4,0,0,0,1,0,0,0,1,0,0,0,1,1,1,4,4,4,4,1,1,1,
     :     1,1,1/
      DATA Z2 /1,2,2,0,1,1,1,1,2,2,2,2, 2,6,6,8, 0, 0,0,0,1,3,4,4,6,7,8,
     :     6,7,8,3,5,6,6,8,8,8,8,9,9,9,10,10,10,10,10,10,11,11,11,11,11,
     :     11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,
     :     13,13,13,13,13,13,13,13,13,13,13,11,11,11,14,14,14,7,8,10/
      CXA = 5.0/3.0
      CXC = CXA - 1.0
      CXB = 2.0*CXC
      CXD = 1.86
      DO J = 1, 92
         CZA(J) = (Z1(J)+Z2(J))**CXA - Z1(J)**CXA - Z2(J)**CXA
         CZB(J) = (Z1(J)+Z2(J))**CXB - Z1(J)**CXB - Z2(J)**CXB
         CZC(J) = -((Z1(J)+Z2(J))**CXC - Z1(J)**CXC - Z2(J)**CXC)
         CZD(J) = (Z1(J)+Z2(J))**CXD - Z1(J)**CXD - Z2(J)**CXD
      END DO
* atomic data
      DATA IZ /1, 2,  6,  7,  8, 10, 12, 14, 26, 2/
      DATA BN /1, 4, 12, 14, 16, 20, 24, 28, 56, 3/
      DATA AM /1.0078, 4.0026, 12.0, 14.003, 15.995, 19.992, 23.985,
     &     27.977, 55.935, 3.0160/
      DATA OMG /1,2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,10,21,28,25,6,25,
     &     28,21/
      DATA CHI/13.595d0, 25*0.0d0,
     &     24.580d0, 54.403d0, 24*0.0d0,
     &     11.26d0,24.38d0,47.86d0,64.48d0,391.99d0,489.84d0,20*0.0d0,
     &     14.54d0,29.60d0,47.43d0,77.45d0,97.86d0,551.92d0,666.83d0,19*0.0d0,
     &     13.61d0,35.15d0,54.93d0,77.39d0,113.87d0,138.08d0,739.11d0,871.12d0,18*0.0d0,
     &     21.56d0,41.07d0,63.5d0,97.16d0,126.4d0,157.91d0,207.3d0,239.0d0,1196.0d0,1360.0d0,
     &     16*0.0d0,
     &     7.64d0,15.03d0,80.12d0,109.29d0,141.23d0,186.86d0,225.31d0,265.96d0,327.90d0,
     &     367.36d0,1761.2d0,2085.0d0,14*0.0d0,     
     &     8.15d0,16.34d0,33.46d0,45.13d0,166.73d0,205.11d0,246.41d0,303.87d0,
     &     351.83d0,401.3d0,476.0d0,523.2d0,2436.0d0,2666.0d0,12*0.0d0,
     &     7.90d0,16.18d0,30.64d0,56.0d0,79.0d0,105.0d0,133.0d0,
     &     151.0d0,235.0d0,262.0d0,290.0d0,321.0d0,
     &     355.0d0,390.0d0,457.0d0,489.0d0,1266.0d0,1358.0d0,1456.0d0,
     &     1582.0d0,1689.0d0,1799.0d0,
     &     1950.0d0,2045.0d0,8828.0d0,9278.0d0/
      DATA DH2/4.477, -0.448, 0.1562, -0.0851/
* atomic masses consistent with Q-values
      CQ = CMEVMU/CL**2
      AM(3) = 12.0
      AM(2) = (AM(3) + CQ*QRT(9))/3.0
      AM(1) = (AM(2) + CQ*(2*QRT(1) + QRT(2)))/4.0
      AM(4) = AM(3) + 2*AM(1) - CQ*QRT(6)
      AM(5) = AM(3) + AM(2) - CQ*QRT(10)
      AM(6) = AM(5) + AM(2) - CQ*QRT(12)
      AM(7) = AM(6) + AM(2) - CQ*QRT(13)
      AM(10) = (AM(2) +2*AM(1) +CQ*QRT(2))/2.0
* sum ionization potentials
C Need to find CHI's for He3
      DO I = 1, 9
         DO J = 2, IZ(I)
            CHI(J,I) = CHI(J,I) + CHI(J-1,I)
         END DO
      END DO
* constants for electron screening
      DO I = 1,10
         VZ(I) = IZ(I)**(3*CXD - 4.0)
      END DO

* Table values for mass loss

      AIJ(1,1)=6.34916
      AIJ(1,2)=-5.04240
      AIJ(1,3)=-0.83426
      AIJ(1,4)=-1.13925
      AIJ(1,5)=-0.12202
      AIJ(2,1)=3.41678
      AIJ(2,2)=0.15629
      AIJ(2,3)=2.96244
      AIJ(2,4)=0.33659
      AIJ(2,5)=0.57576
      AIJ(3,1)=-1.08683
      AIJ(3,2)=0.41952
      AIJ(3,3)=-1.37272
      AIJ(3,4)=-1.07493
      AIJ(4,1)=0.13095
      AIJ(4,2)=-0.09825
      AIJ(4,3)=0.13025
      AIJ(5,1)=0.22427
      AIJ(5,2)=0.46591
      AIJ(6,1)=0.11968


      END
