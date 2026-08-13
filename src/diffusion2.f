C     Compute diffusion coefficients from Paquette et al. 1986
C     For minor species
      SUBROUTINE DIFFUSION2(RHO, T)

      IMPLICIT NONE

      REAL*8 CG, DNSUM, MIZZ, MT, RSUN, CNSTS, ABUND, WW2
      REAL*8 ECHAR, DIFFU2, XA2, PHITT11, MSUM, GAMMATT, SUMNZ2, GAMMA1T
      REAL*8 DC, AVM, E, DH2, EPS1T, F1T22, EPSTT, P1T
      REAL*8 DELTA, PSI1T, CKBA, Q1, DN1, CEN, CEVB, QT
      REAL*8 FTT22, AMU, CLN10, OMG, DLOG, X, CB, Q1T
      REAL*8 LAMBDA, CCHI, PT, PSIPSIN, DNT, M1, C, X1
      REAL*8 XT, CR, PI4, BARYN, PHITT22, S1, LSUN, A
      REAL*8 RHO, ST, BN, MIN, A12, ABUND2, LE, T
      REAL*8 PSIN1PSI, MAX, EQ, CD, XW, F1T1, LAMBDAD, PHI1T22
      REAL*8 LAMBDAI, PSIN, PHI1T1, CC, DMIN1, CMEVMU, CSECYR, DD
      REAL*8 PSIN1, PHI1122, CA, E11, D, FTT11, DIFCOE, CPL
      REAL*8 TSUNYR, ATDATA, AMASS, DMAX1, CPI, MSUN, P1, PSITT
      REAL*8 AUXIN, DEXP, B
      INTEGER INT, JW, IZZ, IML, N, ICN, NE, IOP
      INTEGER INUC, IBC, J, ION, I, ICL

      COMMON /AUXIN / ICL, ION, JW, IOP, INUC, IBC, ICN, IML(2)
      COMMON /ATDATA/ DH2(4), CCHI(26,9), OMG(27), AMASS(10), BN(10), IZZ(10)
      COMMON /ABUND / X(10),XW(13),AVM
      COMMON /ABUND2/ XA2(50), WW2(50)
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CC, CD, CG, CKBA, CR(1), CEVB, 
     &                CEN, CPL, CMEVMU, CSECYR, LSUN, MSUN, RSUN, TSUNYR
      COMMON /DIFFU2/ D(50), A12(50)
      COMMON /DIFCOE/ DC(50,4,3), DD(50,4)
      DIMENSION F1t1(3), PHI1t1(3)
      DIMENSION BARYN(50), MIZZ(50)
      DATA BARYN/61.0, 1.0, 2.0, 3.0, 7.0, 7.0, 11.0, 13.0, 14.0,  15.0,
     : 17.0, 18.0, 19.0, 21.0, 22.0, 22.0, 23.0, 24.0, 25.0, 26.0,
     : 26.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,56.0,57.0,58.0,
     : 59.0,60.0,59.0,58.0,59.0,60.0,61.0,1.0,4.0,12.0,14.0,16.0,20.0,1.0,1.0,1.0,1.0/ 
      DATA MIZZ/30,0,1,2,3,4,5,6,6,7,8,8,9,10,10,11,11,12,12,12,13,13,13,14,14,14,15,
     :16,16,16,26,26,26,26,26,27,28,28,28,28,1,2,6,7,8,10,0,0,0,0/
C Important consts not passed in CNSTS
      AMU = 1.6605402D-24
      ECHAR = 4.8032068D-10
      LAMBDAI = (3d0/(PI4*RHO/(AVM*AMU)))**(1d0/3d0) ! should be n_i...
C Compute Debye length
      LAMBDAD = CKBA*AMU*T/(PI4*(ECHAR**2d0))
      SUMNZ2 = 0d0
      DO I = 1,50
            IF (I.NE.2) THEN
                  SUMNZ2 = SUMNZ2+(MIZZ(I)**2d0)*RHO*XA2(I)/(BARYN(I)*AMU)
            END IF
      END DO

      LAMBDAD = (LAMBDAD/SUMNZ2)**(1d0/2d0)
C pick max of lambdaI and lambdad
      LAMBDA = DMAX1(LAMBDAI,LAMBDAD)

      DO I=1,50
            IF (I.NE.2) THEN
C assume species 1 is hydrogen
                  gamma1t = 4d0*CKBA*AMU*T*LAMBDA/(1d0*MIZZ(I)*ECHAR**2d0)
                  psi1t = DLOG(DLOG(1+gamma1t**(2d0)))
                  msum = 1d0 + BARYN(I)
                  m1 = 1d0/msum
                  mt = BARYN(I)/msum
C     number densities of species
                  dn1 = RHO*XA2(41)/(AMU*BARYN(41))
                  dnt = RHO*XA2(I)/(AMU*BARYN(I))
                  dnsum = dn1 + dnt
                  x1 = dn1/dnsum
                  xt = dnt/dnsum
                  eps1t = CPI*(1d0*MIZZ(I)*ECHAR**2d0/(2d0*CKBA*AMU*T))**2d0
                  eps1t = eps1t*(CKBA*T/(2d0*CPI*(msum*m1*mt)))**(0.5d0)
C Choose appropriate psi region
C note we're only using the repulsive potentials
                  IF (psi1t.LE.3.0) THEN
C Determine appropriate n value for psi1t
                        N = MAX(1, MIN(50, INT((psi1t +7d0)/0.2d0)))
                        PSIN1 = -7d0 + 0.2d0*(N+1)
                        PSIN = -7d0 + 0.2d0*N
                        psin1psi = PSIN1 - psi1t
                        psipsin = psi1t - PSIN
C use splines to calculate collision integrals
                        DO J=1,3
                              F1t1(J) = DEXP(DC(N,1,J)*psin1psi**3.0 + DC(N,2,J)*psipsin**3.0
     +                                + DC(N,3,J)*psin1psi + DC(N,4,J)*psipsin)
                              PHI1t1(J) = F1t1(J)*eps1t
                        END DO

                        F1t22 = DEXP(DD(N,1)*psin1psi**3.0 + DD(N,2)*psipsin**3.0
     +                        + DD(N,3)*psin1psi + DD(N,4)*psipsin)
                        PHI1t22 = F1t22*eps1t
                  ELSE
C     Doesn't need splitting in two as we only consider repulsive potentials
                        F1t1(1) = 1.00141d0*DEXP(psi1t) - 3.18209d0
                        F1t1(2) = 0.99559d0*DEXP(psi1t) - 1.29553d0
                        F1t1(3) = 1.99814d0*DEXP(psi1t) - 0.64413d0
                        F1t22 = 1.99016d0*DEXP(psi1t) - 4.56958d0
                        PHI1t1(1) = F1t1(1)*eps1t
                        PHI1t1(2) = F1t1(2)*eps1t
                        PHI1t1(3) = F1t1(3)*eps1t
                        PHI1t22 = F1t22*eps1t
                  END IF
            
                  IF (I.EQ.41) THEN
C Self-interaction for H
                        PHI1122 = PHI1t22
                        E11 = CKBA*AMU*T/(8d0*m1*m1*PHItt11)
                        P1 = 8d0*m1*E11*PHI1122/(5d0*CKBA*AMU*T)
                  ELSE
C also need self interaction terms for each species
                        gammatt = 4d0*CKBA*AMU*T*LAMBDA/((MIZZ(I)*ECHAR)**2d0)
                        psitt = DLOG(DLOG(1+gammatt**(2d0)))
                        epstt = CPI*((MIZZ(I)*ECHAR)**2d0/(2d0*CKBA*AMU*T))**2d0
                        epstt = epstt*(CKBA*T/(2d0*CPI*(0.5*BARYN(I))))**(0.5d0)
                        IF (psitt.LE.3.0) THEN
C Determine appropriate n value for psi1t
                              N = MAX(1, MIN(50, INT((psitt +7d0)/0.2d0)))
                              PSIN1 = -7d0 + 0.2d0*(N+1)
                              PSIN = -7d0 + 0.2d0*N
                              psin1psi = PSIN1 - psitt
                              psipsin = psitt - PSIN
                              Ftt11 = DEXP(DC(N,1,1)*psin1psi**3.0 + DC(N,2,1)*psipsin**3.0
     +                              + DC(N,3,1)*psin1psi + DC(N,4,1)*psipsin)
                              Ftt22 = DEXP(DD(N,1)*psin1psi**3.0 + DD(N,2)*psipsin**3.0
     +                              + DD(N,3)*psin1psi + DD(N,4)*psipsin)
                              PHItt11 = Ftt11*epstt
                              PHItt22 = Ftt22*epstt
                        ELSE
                              Ftt22 = 1.99016d0*DEXP(psitt) - 4.56958d0
                              PHItt22 = Ftt22*epstt
                        END IF
                  END IF
C         if (I.EQ.2) write (*,*) I, PHI1t1(1), PHI1t1(2), PHI1t1(3), PHI1t22
C         if (I.EQ.2) write(*,*) psitt, PHItt22, epstt, eps1t
                  A = PHI1t22/(5d0*PHI1t1(1))
                  B = (5d0*PHI1t1(2) - PHI1t1(3))/(5d0*PHI1t1(1))
                  C = (2d0*PHI1t1(2)/(5d0*PHI1t1(1))) - 1d0
C     if (I.EQ.2) write (*,*) I, A, B, C
                  E = CKBA*AMU*T/(8d0*m1*mt*PHI1t1(1))
                  P1t = 3d0*(m1 - mt)**2d0 +4d0*m1*mt*A
C Is E right for P_i? It's not clear what subscripts should be on omega, 
C which for some reason I've called PHI. Don't know what I was on...
C original implementation
                  P1 = 8d0*m1*E*PHI1122/(5d0*CKBA*AMU*T)
                  Pt = 8d0*mt*E*PHItt22/(5d0*CKBA*AMU*T)
C New implementation -- assumes E should be the self-property
C         Ett = CKBA*AMU*T/(8d0*mt*mt*PHItt11)
C         Pt = 8d0*mt*Ett*PHItt22/(5d0*CKBA*AMU*T)
C         if (I.EQ.2) write (*,*) P1t, P1, Pt
                  Q1 = P1*(6d0*mt**2d0 + 5*m1**2d0 - 4d0*(m1**2d0)*B + 8*m1*mt*A)
                  Qt = Pt*(6d0*m1**2d0 + 5*mt**2d0 - 4d0*(mt**2d0)*B + 8*mt*m1*A)
                  Q1t = 3d0*(m1 - mt)**2d0
                  Q1t = Q1t*(5d0-4d0*B) + 4*m1*mt*A*(11d0-4d0*B)+2d0*P1*Pt
C         if (I.EQ.2) write (*,*) Q1t, Q1, Qt
                  S1 = m1*P1 - mt*(3d0*(mt-m1)+4d0*m1*A)
                  St = mt*Pt - m1*(3d0*(m1-mt)+4d0*mt*A)
                  delta = 5d0*(C**2d0)*((m1**2d0)*P1*(x1**2d0) + (mt**2d0)*Pt*(xt**2d0)
     :           + P1t*x1*xt)
                  delta = delta/((x1**2d0)*Q1 + (xt**2d0)*Qt + x1*xt*Q1t)
C delta fudge to stop getting silly values
                  delta = dmin1(delta, 0.1d0)
C         if (I.EQ.2) write (*,*) S1, St, delta
C         IF (delta.gt.1d0) write (*,*) "delta problems!", delta, I
C         IF (delta.gt.1d0.AND.I.EQ.2) THEN
C                 write (*,*) A, B, T, RHO, X(1), X(2)
C                 write (*,*) ((x1**2d0)*Q1 + (xt**2d0)*Qt + x1*xt*Q1t)
C                 write (*,*) 5d0*C**2d0*((m1**2d0)*P1*(x1**2d0) + (mt**2d0)*Pt*(xt**2d0)
C     :        + P1t*x1*xt)
C                 write (*,*) C, m1, mt, x1, xt, P1, Pt, P1t, Q1, Qt, Q1t
C                 stop
C         END IF
                  D(I) = 3d0*CKBA*T/(16d0*dnsum*msum*m1*mt*PHI1t1(1))
C second approx to diffusion coeff.
                  D(I) = D(I)/(1-delta)
C thermodiffusion coeff. -- if I need it...
                  A12(I) = 5d0*C*(x1*S1 - xt*St)/((x1**2d0)*Q1 + (xt**2d0)*Qt + x1*xt*Q1t)
            END IF
      END DO
      RETURN
      END
