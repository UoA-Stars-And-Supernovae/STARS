      SUBROUTINE EQUNS1(K, K1, K2, ISTAR, IVAR)

      IMPLICIT NONE

      REAL*8 VP, WI, HT, D14, MT, D3, X16T, ETH
      REAL*8 ACCRET, D4, OUTE, VX, DABS, D16, SG2, GTA
      REAL*8 TRC2, MU, MWT, SG1, VT, SGTH, LK, DA16
      REAL*8 SG, EQU, VMK, TRC1, WT2, DA4, X20, DA20
      REAL*8 D12, X12, X24, DA12, D20, GRAD, EGR, X1
      REAL*8 X3, X1T, VTK, QK, GT, ACCOMPOS, X14, VARACC
      REAL*8 EVMODE, X12T, BCHORB, X3T, R2K, VM, X4, LE
      REAL*8 T, X4T, OP, WT3, EQ, BC3, SGLEV, LQ
      REAL*8 R, BC2, BC1, VVM, DT2, L, DA3, MWTS
      REAL*8 WTM, BCHSPIN, DD, RLF, X16, VPK, GR, DT1
      REAL*8 MESH, TANH, ZS, R2, X20T, LT, LEDD, DA14
      REAL*8 TRANS, X14T, QQ, PS
      INTEGER NMAXMSH, K1, IMODE, K2, NE, IOTHER, IVARACC, K
      INTEGER INE, IVAR, ISTAR, IVMC, IVMS, IMLWR

      PARAMETER (NMAXMSH = 2000)
      COMMON /INE   / BC1(3), BC2(3), BC3(3), BCHORB(3), BCHSPIN(3),
     :            VP(3), VPK(3), R2(3), 
     :            R2K(3), VT(3), VTK(3), L(3), LK(3), LQ(3), GTA(3), 
     :            MT(3), VM(3), VMK(3), QK(3), SG(3), T(3), X1(3),  
     :            X1T(3), X16(3), X16T(3), X4(3), X4T(3), X12(3),  
     :            X12T(3), X20(3), X20T(3), X14(3), X14T(3), X24(3),
     :            SGTH(3), MU(3), X3(3), X3T(3), SGLEV(3), DA4(3),
     :            DA12(3), DA14(3), DA16(3), DA20(3), DA3(3),
     :            D4(3), D12(3), D14(3), D16(3), D20(3), D3(3),
     :            WI(309)
      COMMON /OP    / ZS, LEDD, VVM, GR, GRAD, ETH, RLF, EGR, R, QQ
C VAR(3),(2),(1) are values of VAR at current, previous and anteprevious meshpts
      COMMON /OUTE  / EQU(50)
      COMMON /MESH  / TRC1,TRC2,DD,DT1,DT2,MWT,MWTS, IVMC,IVMS
      COMMON /ACCRET/ ACCOMPOS(7,31, 2)
      COMMON /TRANS / HT(26,NMAXMSH,2)
      COMMON /EVMODE/ IMODE
      COMMON /VARACC/ IVARACC, IMLWR

      PS(VX) = 0.5D0*(VX+DABS(VX))
C 30/5/03 RJS Smooth viscous mesh
      WTM = 0.5 + 0.5*TANH((K - TRC1)/1.5)
      WTM = MWT*WTM

C Surface mesh viscosity
      IF (IVMS.EQ.1) THEN
            WTM = WTM + MWTS*(0.5 - 0.5*TANH((K - TRC2)/1.5))
      END IF

      IF (WTM.GT.1.0) THEN
            WTM = 1.0
      ELSE IF (WTM.LT.0.0) THEN
            WTM = 0.0
      END IF

      IF ( K.LE.K1 ) THEN
C surface boundary conditions
            EQU(1) = BC1(3)
            EQU(2) = BC2(3)
            EQU(3) = BC3(3)
C Orbital angular momentum
            EQU(4) = BCHORB(3)
C Spin period of star
            EQU(5) = BCHSPIN(3)

            RETURN
      ELSE IF ( K.LE.K2 ) THEN
C first-order difference equations at interior points
            WT3 = 0.5D0
C weighted alternative to central differencing
C        WT3 = 0.5D0*WT(3)/(1.0D0+WT(3))
            WT2 = 1.0D0 - WT3
            EQU(1) = VP(3) - VP(2) - WT3*VPK(3)-WT2*VPK(2)
            EQU(2) = R2(3) - R2(2) - WT2*R2K(3)-WT3*R2K(2)
            EQU(3) = VT(3) - VT(2) - WT3*VTK(3)-WT2*VTK(2)
            EQU(4) = L(3) - L(2) - WT2*LK(3)-WT3*LK(2)
     :              - LQ(2)*GTA(2)*PS(MT(2)) + LQ(3)*GTA(3)*PS(-MT(3))
            EQU(5) = VM(3) - VM(2) - 0.5D0*(VMK(3)+VMK(2))
C 22/3/03 RJS Added viscous mesh
            EQU(6)=(1.0-WTM)*(QK(3) - QK(2))+3.0d7*WTM*MT(3)
            EQU(13) = BCHSPIN(3) - BCHSPIN(2)

            IF ( K.EQ.K1+1 ) THEN
C next-to-surface boundary conditions for second-order equations
C Attempt at variable composition accretion - only if in binary mode
                  IOTHER = 3 - ISTAR ! This maps 1 -> 2 and 2 -> 3
                  IF ((HT(24,1,IOTHER).GT.0d0.OR.HT(23,1,ISTAR).LT.0d0).AND.IMODE.EQ.2 .AND. IVARACC.EQ.1) THEN
C              If both stars are filling their roche lobes, set the abundance of the accreted material
C              to the average of the two stars
C                IF (HT(24,1,IOTHER).GT.0d0.AND.HT(24,1,ISTAR).GT.0d0.AND.IVARACC.EQ.1) THEN
                        EQU(7) = 0.5*(ACCOMPOS(1,IVAR+1,IOTHER)+ X1(3)) - X1(2)
                        EQU(8) = 0.5*(ACCOMPOS(5,IVAR+1,IOTHER) + X16(3)) - X16(2)
                        EQU(9) = 0.5*(ACCOMPOS(2,IVAR+1,IOTHER) + X4(3)) - X4(2)
                        EQU(10) = 0.5*(ACCOMPOS(3,IVAR+1,IOTHER) + X12(3)) - X12(2)
                        EQU(11) = 0.5*(ACCOMPOS(6,IVAR+1,IOTHER) + X20(3)) - X20(2)
                        EQU(12) = 0.5*(ACCOMPOS(4,IVAR+1,IOTHER) + X14(3)) - X14(2)
                        EQU(14) = 0.5*(ACCOMPOS(7,IVAR+1,IOTHER) + X3(3)) - X3(2)
                ! SMR + JJE 21 November 2023
C The following lines are a more accurate but less stable implementation of the
C next to surface boundary condition -- should be used if you want to do thermohaline mixing
C                SG2 = -(PS(MT(2))+1d-5) !1d-5
C                EQU(7) = SG2*(X1(3)-X1(2)) + PS(MT(2))*(X1(2)-ACCOMPOS(1,IVAR+1,IOTHER))
C      :              - X1T(2)
C                EQU(8) = SG2*(X16(3)-X16(2)) + PS(MT(2))*(X16(2)-ACCOMPOS(5,IVAR+1,IOTHER))
C      :              - X16T(2)
C                EQU(9) = SG2*(X4(3)-X4(2)) + PS(MT(2))*(X4(2)-ACCOMPOS(2,IVAR+1,IOTHER))
C      :              - X4T(2)
C                EQU(10) = SG2*(X12(3)-X12(2)) + PS(MT(2))*(X12(2)-ACCOMPOS(3,IVAR+1,IOTHER))
C      :              - X12T(2)
C                EQU(11) = SG2*(X20(3)-X20(2)) + PS(MT(2))*(X20(2)-ACCOMPOS(6,IVAR+1,IOTHER))
C      :              - X20T(2)
C                EQU(12) = SG2*(X14(3)-X14(2)) + PS(MT(2))*(X14(2)-ACCOMPOS(4,IVAR+1,IOTHER))
C      :              - X14T(2)
C                EQU(14) = SG2*(X3(3)-X3(2)) + PS(MT(2))*(X3(2)-ACCOMPOS(7,IVAR+1,IOTHER))
C      :              - X3T(2)
                  ELSE
                        EQU(7) = X1(3) - X1(2)
                        EQU(8) = X16(3) - X16(2)
                        EQU(9) = X4(3) - X4(2)
                        EQU(10) = X12(3) - X12(2)
                        EQU(11) = X20(3) - X20(2)
                        EQU(12) = X14(3) - X14(2)
                        EQU(14) = X3(3) - X3(2)
                  END IF

                  RETURN
            ELSE
C second-order difference equations at interior points
                  SG1 = 0.5D0*(SG(1)+SG(2)) - PS(MT(2))
                  SG2 = 0.5D0*(SG(2)+SG(3)) - PS(-MT(3))
C Add in thermohaline mixing
                  SG1 = SG1 + 0.5*(SGTH(1)+SGTH(2))*PS((MU(1)-MU(2)))
                  SG2 = SG2 + 0.5*(SGTH(2)+SGTH(3))*PS((MU(2)-MU(3)))
C Note gravitational settling is hard-wired into having H as the dominant element!
                  EQU(7) = (SG2 + 0.5*(DA4(2)+DA4(3)))*(X1(3)-X1(2))
     :                   - (SG1+0.5*(DA4(1)+DA4(2)))*(X1(2)-X1(1)) - X1T(2)
     :                   + D4(2)*X4(2) - D4(3)*X4(3)
     :                   + D12(2)*X12(2) - D12(3)*X12(3)
     :                   + D14(2)*X14(2) - D14(3)*X14(3)
     :                   + D16(2)*X16(2) - D16(3)*X16(3)
     :                   + D20(2)*X20(2) - D20(3)*X20(3)
     :                   + D3(2)*X3(2) - D3(3)*X3(3)
                  EQU(8) = (SG2 + 0.5*(DA16(2)+DA16(3)))*(X16(3)-X16(2))
     :                   - (SG1+0.5*(DA16(1)+DA16(2)))*(X16(2)-X16(1)) - X16T(2)
     :                   - D16(2)*X16(2) + D16(3)*X16(3)
                  EQU(9) = (SG2 + 0.5*(DA4(2)+DA4(3)))*(X4(3)-X4(2))
     :                   - (SG1+0.5*(DA4(1)+DA4(2)))*(X4(2)-X4(1)) - X4T(2)
     :                   - D4(2)*X4(2) + D4(3)*X4(3)
                  EQU(10) = (SG2 + 0.5*(DA12(2)+DA12(3)))*(X12(3)-X12(2))
     :                   - (SG1+0.5*(DA12(1)+DA12(2)))*(X12(2)-X12(1)) - X12T(2)
     :                   - D12(2)*X12(2) + D12(3)*X12(3)
                  EQU(11)= (SG2 + 0.5*(DA20(2)+DA20(3)))*(X20(3)-X20(2))
     :                   - (SG1+0.5*(DA20(1)+DA20(2)))*(X20(2)-X20(1)) - X20T(2)
     :                   - D20(2)*X20(2) + D20(3)*X20(3)
                  EQU(12)= (SG2 + 0.5*(DA14(2)+DA14(3)))*(X14(3)-X14(2))
     :                   - (SG1+0.5*(DA14(1)+DA14(2)))*(X14(2)-X14(1)) - X14T(2)
     :                   - D14(2)*X14(2) + D14(3)*X14(3)
                  EQU(14)= (SG2 + 0.5*(DA3(2)+DA3(3)))*(X3(3)-X3(2))
     :                   - (SG1+0.5*(DA3(1)+DA3(2)))*(X3(2)-X3(1)) - X3T(2)
     :                   - D3(2)*X3(2) + D3(3)*X3(3)
                  RETURN
            END IF
      END IF
C central boundary conditions for first-order equations
      IF (WTM.GT.1.0) THEN
            WTM = 1.0
      END IF

      IF (WTM.LT.0.0) THEN
            WTM = 0.0
      END IF

      IF (WTM.NE.0.0) THEN
            EQU(1) = MT(3)
      ELSE
            EQU(1) = VM(3) + 1.5D0*VMK(3) - 0.5D0*VMK(2)
      END IF
C This was the original central L boundary condition -- PPE said there
C was a good reason for it but he couldn't remember it.
C      EQU(2) = L(3) + 0.93333D0*LK(3) - 0.1885D0*LK(2) - LQ(3)*GTA(3)
C     :         *MT(3)
      EQU(2) = L(3) + 1.5D0*LK(3) - 0.5D0*LK(2) - LQ(3)*GTA(3)*MT(3)
      EQU(3) = R2(3) + 1.5D0*R2K(3) - 0.5D0*R2K(2)
C central boundary conditions for second-order equations
      SG2 = 0.5D0*(SG(2)+SG(3))
C Plus thermohaline mixing
      SG2 = SG2 + 0.5*(SGTH(2)+SGTH(3))*PS((MU(2)-MU(3)))
      EQU(4) = SG2*(X1(3)-X1(2)) + X1T(3)
      EQU(5) = SG2*(X16(3)-X16(2)) + X16T(3)
      EQU(6) = SG2*(X4(3)-X4(2)) + X4T(3)
      EQU(7) = SG2*(X12(3)-X12(2)) + X12T(3)
      EQU(8) = SG2*(X20(3)-X20(2)) + X20T(3)
      EQU(9) = SG2*(X14(3)-X14(2)) + X14T(3)
      EQU(10) = SG2*(X3(3)-X3(2)) + X3T(3)

      RETURN
      END
