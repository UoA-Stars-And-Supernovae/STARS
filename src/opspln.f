      SUBROUTINE OPSPLN

      IMPLICIT NONE

      REAL*8 CSX, CS, DFLOAT, W, STAT1, STAT3, FSPL, FKLM
      REAL*8 FKHM, TFM, FRM, CNIU
      INTEGER MT, JR, IQ, IT, MR, IC, JT
      INTEGER JQ, IR, IFAIL, JCSX, JX

      PARAMETER (MT = 127,MR = 90)
      REAL*8 MAT(4,MT)
      COMMON /STAT1 / CSX(10), CS(90,127,10), CNIU(60,41,2), W(18400), 
     :                JCSX
      COMMON /STAT3/ FSPL(4,4,127,90,10),TFM(127),FRM(90),
     :                FKLM(6),FKHM(6)
*
* Calculate a bicubic spline interpolation fit for the temperature
* and density opacity fit.  This is based on routines from Princeton
* Splines, D. McCune found by L. Dray.
*
      DO 10 JT = 1,MT
      TFM(JT) = 2.95D0 + 0.05D0*DFLOAT(JT)
   10 CONTINUE
      DO 20 JR = 1,MR
      FRM(JR) = 0.25D0*DFLOAT(JR) - 12.25D0
   20 CONTINUE
      IFAIL = 0
      DO 30 JX = 1,10
      DO 26 JQ = 1,MT
      DO 23 IQ = 1,MR
      FSPL(1,1,JQ,IQ,JX) = CS(IQ,JQ,JX)
   23 CONTINUE        
   26 CONTINUE
*
* Construct splines in the T direction.
*
      DO IR=1,MR
         DO IT = 1,MT
            MAT(1,IT) = FSPL(1,1,IT,IR,JX)
         ENDDO
         CALL SPLINE(MT,TFM,MAT)
         DO IT = 1,MT-1
            FSPL(2,1,IT,IR,JX) = MAT(2,IT)
            FSPL(3,1,IT,IR,JX) = MAT(3,IT)
            FSPL(4,1,IT,IR,JX) = MAT(4,IT)
         ENDDO
*
      ENDDO
*
* Construct splines in the rho direction.
*
      DO IT = 1,MT-1
*
* Construct splines for each T coeff
*
         DO IC = 1,4
            DO IR = 1,MR
               MAT(1,IR) = FSPL(IC,1,IT,IR,JX)
            ENDDO
            MAT(2,1) = 0D0
            MAT(3,1) = 0D0
            MAT(2,MR)= 0D0
            MAT(3,MR)= 0D0
            CALL SPLINE(MR,FRM,MAT)
            DO IR = 1,MR-1
               FSPL(IC,2,IT,IR,JX) = MAT(2,IR)
               FSPL(IC,3,IT,IR,JX) = MAT(3,IR)
               FSPL(IC,4,IT,IR,JX) = MAT(4,IR)
            ENDDO
C
         ENDDO
C
      ENDDO
C
   30 CONTINUE
      RETURN
      END
