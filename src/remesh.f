**==REMESH.FOR
      SUBROUTINE REMESH(M, NCH, CH, CO, CC, CNE )
      IMPLICIT NONE
      REAL*8 CG, XH, BM, XO17, FAM, RSUN, OS, XS32
      REAL*8 CNSTS, XAL27, AM, XS33, XCO59, XMG25, VX, TRB
      REAL*8 RESETTING, XN15, XBE, XNE, CHE, DK, XB11, EP
      REAL*8 XSI, XNE21, XNE22, XO, AGE, DH2, X26M, Q2
      REAL*8 XN, HNUC, XNA23, XSI29, Q1, XSI30, CEN, Q
      REAL*8 CEVB, RMG, CLN10, XNI58, XC13, SPIN, TC, OMG
      REAL*8 XC14, DLOG, CH, DVAR, CB, CO, DT, XHE3
      REAL*8 AMTB, AMTA, XMG, AM0, VROT1, VAR, XMG26, RHL
      REAL*8 DH, CT, XNA22, XFE60, XP31, OMEGA, CNE, AK2
      REAL*8 CR, PI4, XFE56, XNI61, CSI, XFE59, CLSUN, CLIGHT
      REAL*8 XS34, XFE57, DR, XNI60, W, XHE, PHI, EVMODE
      REAL*8 XFE, YUK, CN, ORBITAL, VMA, BN, TM, VI
      REAL*8 DINF, XLI, QE, XF, OP, ALPHA, EQ, CD
      REAL*8 CM, X26G, RML, T0, AM1, XO18, REMESHING, VROT
      REAL*8 CC, DHNUC, CFE, CMEVMU, CSECYR, ANG, CMG, CA
      REAL*8 XSI28, RCD, XMG24, CPL, GE, H, FMAC, TSUNYR
      REAL*8 ATDATA, XF19, ANGMOM, CPI, EC, CHI, XD, VROT2
      REAL*8 CMSUN, XNI59, ZS, RMT, XC, XFE58, DSQRT, SODDS
      REAL*8 QQ, DEXP, AK1, AC
      INTEGER JIN, IZ, KK, IRS1, IMODE, M, NMESHORIG, ISX
      INTEGER MODEL, IB, NE, IRAM, NMESH, K, NCH, ID
      INTEGER IRS2, JZ, NUCMAT, MAXMSH, J, I, INERTI, ISTAR
      INTEGER INF, LT

      PARAMETER (MAXMSH = 2000)
      COMMON H(60,MAXMSH),DH(60,MAXMSH),EP(3),NMESH,JIN,ID(200)
      COMMON /NUCMAT/ HNUC(100,MAXMSH), DHNUC(100,MAXMSH)
      COMMON /SODDS / ALPHA, RML, CMG, CSI, CFE, CT(10), AGE, DT, AM1,
     :  EC, BM, ANG, CM, AMTA, AMTB, TM(2), T0, AM0, TC(2), OS, AC, RCD,  
     :  RMG, RHL, XF, DR, AK1 ,RMT, AK2, IZ(4), IB, ISX(45),
     :  TRB
      COMMON /INF   / VAR(60)
      COMMON /DINF  / DVAR(60)
      COMMON /OP    / ZS, W(8), QQ
      COMMON /ATDATA/ DH2(4), CHI(26,9), OMG(27), AM(10), BN(10), JZ(10)
      COMMON /YUK   / Q(MAXMSH)
      COMMON /EVMODE/ IMODE
      COMMON /CNSTS / CPI, PI4, CLN10, CA, CB, CLIGHT, CD, CG, CR(2), CEVB, 
     &                CEN, CPL, CMEVMU, CSECYR, CLSUN, CMSUN, RSUN, TSUNYR
      COMMON /INERTI/ VI(2)
      COMMON /ANGMOM/ VROT1, VROT2, FMAC, FAM, IRAM, IRS1, IRS2
      WRITE(32,*) "Remeshing model"
C RJS 22/11/07 - Save NMESH, we need it for the second remesh pass
      NMESHORIG = NMESH
      DO ISTAR = 1,IMODE!2
         NMESH = NMESHORIG
C Calculate moment of inertia
         VI(ISTAR) = 0d0
         DO K=1,NMESH - 1
            VI(ISTAR) = VI(ISTAR) + 2.0/3.0*(DEXP(H(4+15*(ISTAR-1),K)) - 
     :           DEXP(H(4+15*(ISTAR-1),K+1)))*DEXP(H(7+15*(ISTAR-1),K))**2.0
         END DO
         VI(ISTAR) = VI(ISTAR) + 2.0/5.0*DEXP(H(4+15*(ISTAR-1),NMESH))*DEXP(H(7+15*(ISTAR-1),NMESH))**2.0
         IF ( NCH.LT.3 ) GO TO 41
*
* Set initial composition.
* The composition variables are NOT the actual mass fractions if we use
* non-integer atomic masses, so we have to compute what they are
*
         CHE = 1.0 - CH - ZS
         CN = 1D0 - CC - CO - CNE - CMG - CSI - CFE
         XH = CH*BN(1)/AM(1)
         XHE = CHE*BN(2)/AM(2)
         XC = CC*ZS*BN(3)/AM(3)
         XN = CN*ZS*BN(4)/AM(4)
         XO = CO*ZS*BN(5)/AM(5)
         XNE = CNE*ZS*BN(6)/AM(6)
C should do He3 like this...
         XMG = CMG*ZS
         XSI = CSI*ZS
         XFE = CFE*ZS
         VMA = XH + XHE + XC + XN + XO + XNE + XMG + XSI + XFE
C Minor element initial abundances from Anders & Grevesse 89
         XD = 4.801d-5*(ZS/0.02)
         XHE3 = 2.929d-5*(ZS/0.02)
         XLI = 9.353d-9*(ZS/0.02)
C No data available for Be-7, assume it's zero?
         XBE = 0d0
         XB11 = 4.725d-9*(ZS/0.02)
         XC13 = 3.650d-5*(ZS/0.02)
         XC14 = 0d0
         XN15 = 4.363d-6*(ZS/0.02)
         XO17 = 3.887d-6*(ZS/0.02)
         XO18 = 2.167d-5*(ZS/0.02)
         XF19 = 4.051d-7*(ZS/0.02)
         XNE21 = 4.127d-6*(ZS/0.02)
         XNE22 = 1.302d-4*(ZS/0.02)
         XNA22 = 0d0
         XNA23 = 3.339d-5*(ZS/0.02)
         XMG24 = 5.148d-4*(ZS/0.02)
         XMG25 = 6.766d-5*(ZS/0.02)
         XMG26 = 7.760d-5*(ZS/0.02)
         X26M = 0d0
         X26G = 0d0
         XAL27 = 5.798d-5*(ZS/0.02)
         XSI28 = 6.530d-4*(ZS/0.02)
         XSI29 = 3.426d-5*(ZS/0.02)
         XSI30 = 2.352d-5*(ZS/0.02)
         XP31 = 8.155d-5*(ZS/0.02)
         XS32 = 3.958d-4*(ZS/0.02)
         XS33 = 3.222d-6*(ZS/0.02)
         XS34 = 1.866d-5*(ZS/0.02)
         XFE56 = 1.169d-3*(ZS/0.02)
         XFE57 = 2.855d-5*(ZS/0.02)
         XFE58 = 3.697d-6*(ZS/0.02)
         XFE59 = 0d0
         XFE60 = 0d0
         XCO59 = 3.358d-6*(ZS/0.02)
         XNI58 = 4.944d-5*(ZS/0.02)
         XNI59 = 0d0
         XNI60 = 1.958d-5*(ZS/0.02)
         XNI61 = 8.594d-7*(ZS/0.02)
         DO 4 K = 1, NMESH
            H(5+15*(ISTAR-1), K) = XH/VMA
            H(3+15*(ISTAR-1), K) = XO/VMA
            H(9+15*(ISTAR-1), K) = XHE/VMA
            H(10+15*(ISTAR-1),K) = XC/VMA
            H(11+15*(ISTAR-1),K) = XNE/VMA
            H(12+15*(ISTAR-1),K) = XN/VMA
C HE3 for structure code
            H(15+15*(ISTAR - 1),K) = XHE3/VMA
C Beginning of minor variables, in mass order (except gallinoes)
            HNUC(1+50*(ISTAR-1),K) = 0d0 ! gallinoes
            HNUC(2+50*(ISTAR-1),K) = 0d0 ! neutrons
            HNUC(3+50*(ISTAR-1),K) = XD/VMA
            HNUC(4+50*(ISTAR-1),K) = XHE3/VMA
            HNUC(5+50*(ISTAR-1),K) = XLI/VMA
            HNUC(6+50*(ISTAR-1),K) = XBE/VMA
            HNUC(7+50*(ISTAR-1),K) = XB11/VMA
            HNUC(8+50*(ISTAR-1),K) = XC13/VMA
            HNUC(9+50*(ISTAR-1),K) = XC14/VMA
            HNUC(10+50*(ISTAR-1),K) = XN15/VMA
            HNUC(11+50*(ISTAR-1),K) = XO17/VMA
            HNUC(12+50*(ISTAR-1),K) = XO18/VMA
            HNUC(13+50*(ISTAR-1),K) = XF19/VMA
            HNUC(14+50*(ISTAR-1),K) = XNE21/VMA
            HNUC(15+50*(ISTAR-1),K) = XNE22/VMA
            HNUC(16+50*(ISTAR-1),K) = XNA22/VMA
            HNUC(17+50*(ISTAR-1),K) = XNA23/VMA
            HNUC(18+50*(ISTAR-1),K) = XMG24/VMA
            HNUC(19+50*(ISTAR-1),K) = XMG25/VMA
            HNUC(20+50*(ISTAR-1),K) = XMG26/VMA
            HNUC(21+50*(ISTAR-1),K) = X26M/VMA
            HNUC(22+50*(ISTAR-1),K) = X26G/VMA
            HNUC(23+50*(ISTAR-1),K) = XAL27/VMA
            HNUC(24+50*(ISTAR-1),K) = XSI28/VMA
            HNUC(25+50*(ISTAR-1),K) = XSI29/VMA
            HNUC(26+50*(ISTAR-1),K) = XSI30/VMA
            HNUC(27+50*(ISTAR-1),K) = XP31/VMA
            HNUC(28+50*(ISTAR-1),K) = XS32/VMA
            HNUC(29+50*(ISTAR-1),K) = XS33/VMA
            HNUC(30+50*(ISTAR-1),K) = XS34/VMA
            HNUC(31+50*(ISTAR-1),K) = XFE56/VMA
            HNUC(32+50*(ISTAR-1),K) = XFE57/VMA
            HNUC(33+50*(ISTAR-1),K) = XFE58/VMA
            HNUC(34+50*(ISTAR-1),K) = XFE59/VMA
            HNUC(35+50*(ISTAR-1),K) = XFE60/VMA
            HNUC(36+50*(ISTAR-1),K) = XCO59/VMA
            HNUC(37+50*(ISTAR-1),K) = XNI58/VMA
            HNUC(38+50*(ISTAR-1),K) = XNI59/VMA
            HNUC(39+50*(ISTAR-1),K) = XNI60/VMA
            HNUC(40+50*(ISTAR-1),K) = XNI61/VMA !X26T/VMA
            HNUC(41+50*(ISTAR-1),K) = XH/VMA
            HNUC(42+50*(ISTAR-1),K) = XHE/VMA
            HNUC(43+50*(ISTAR-1),K) = XC/VMA
            HNUC(44+50*(ISTAR-1),K) = XN/VMA
            HNUC(45+50*(ISTAR-1),K) = XO/VMA
            HNUC(46+50*(ISTAR-1),K) = XNE/VMA
 4       CONTINUE
 41      CONTINUE
C Orbital stuff doesn't work properly - copy surface point to all others
         DO K=2,NMESH
            H(13,K) = H(13,1)
            H(14,K) = H(14,1)
            H(29,K) = H(29,1)
         END DO
C Reset orbital AM using period in modin
         IF (IRAM.EQ.1.AND.ISTAR.EQ.1) THEN
            WRITE(32,*) "Resetting orbital AM"
            DO K=1,NMESH
               H(13,K) = ANG
            END DO
         END IF
         IF (IRS1.EQ.1.AND.ISTAR.EQ.1) THEN
            WRITE(32,*) "Resetting *1 spin"
            VROT = VROT1*1d3
            OMEGA = VROT/(1d9*DEXP(H(7,1)))
C Convert omega to code units
            OMEGA = OMEGA/DSQRT(CG)
            DO K=1,NMESH
               H(14,K) = VI(1)*OMEGA 
            END DO
         END IF
         IF (IRS2.EQ.1.AND.ISTAR.EQ.2) THEN
            WRITE(32,*) "Resetting *2 spin"
            VROT = VROT2*1d3
            OMEGA = VROT/(1d9*DEXP(H(22,1)))
            OMEGA = OMEGA/DSQRT(CG)
            DO K=1,NMESH
               H(29,K) = VI(2)*OMEGA 
            END DO
         END IF
         QE = 0D0
         DO K = 1, NMESH
            DO J = 1, 15        !60
               DVAR(J) = 0.0D0
               VAR(J) = H(J+15*(ISTAR-1), K)
            END DO
            VAR(6) = 1.0D0
            CALL FUNCS1(-1, K, 1, NMESH, ISTAR)
            Q(K) = QQ
         END DO
         PHI = Q(NMESH) - Q(1)
         Q1 = Q(1)
         Q2 = (M-1.0D0)/PHI
         IF ( NCH.GE.2 ) THEN
C IF REQUIRED, REDISTRIBUTE THE MESH POINTS, EITHER BY CHANGING THEIR NUMBER
C OR BY CHANGING THE MESH-SPACING FUNCTION.
            DO K = 1, NMESH
               VX = H(5+15*(ISTAR-1), K) + 1.0D-10
               H(5+15*(ISTAR-1), K) = DLOG(VX)
               Q(K) = (Q(K)-Q1)*Q2 + 1.0D0
            END DO
            KK = 1
            DO K = 1, M
               DK = 0.0D0
               IF ( K.EQ.M ) KK = NMESH
               IF ( K.NE.1 .AND. K.NE.M ) THEN
                  DO I = 1, 50
                     IF ( K.GE.Q(KK+1) ) KK = KK + 1
                     IF ( K.LT.Q(KK+1) ) GO TO 10
                  END DO
 10               CONTINUE
                  DK = (K-Q(KK))/(Q(KK+1)-Q(KK))
               END IF
               DO J = 1, 15 !60
                  DH(J+15*(ISTAR-1), K) = H(J+15*(ISTAR-1), KK) +
     :                 DK*(H(J+15*(ISTAR-1),KK+1)-H(J+15*(ISTAR-1),KK))
               END DO
C Remesh nucleosynthesis matrix
               DO J=1,50
                  DHNUC(J+50*(ISTAR-1),K) = HNUC(J+50*(ISTAR-1), KK) +
     :                 DK*(HNUC(J+50*(ISTAR-1), KK+1)-HNUC(J+50*(ISTAR-1), KK))
               END DO
            END DO
            NMESH = M
            DO K = 1, NMESH
               DH(5+15*(ISTAR-1), K) = DEXP(DH(5+15*(ISTAR-1),K)) - 1.0D-10
C Is this line still necessary?
               IF ( DH(5+15*(ISTAR-1),K).LT.1.0D-5 ) DH(5+15*(ISTAR-1), K) = 0.0D0
               DO J = 1, 15 !60
                  H(J+15*(ISTAR-1), K) = DH(J+15*(ISTAR-1), K)
                  DH(J+15*(ISTAR-1), K) = 0.0D0
               END DO
               DO J=1,50
                  HNUC(J+50*(ISTAR-1), K) = DHNUC(J+50*(ISTAR-1),K)
                  DHNUC(J+50*(ISTAR-1),K) = 0d0
               END DO
            END DO
         END IF
         DO K = 1, NMESH
            H(6+15*(ISTAR-1), K) = 1.0D0/Q2
         END DO
      END DO
      RETURN
      END
