      SUBROUTINE OPACTY(JX, TF, FR, FKL, FKH)

      IMPLICIT NONE

      REAL*8 CSX, FKLO, FKL, DR, CS, FKHO, W, STAT1
      REAL*8 F, DT, GE, STAT3, TF, FR, FKH, FKLM
      REAL*8 FKHM, TFM, FRM, CNIU
      INTEGER MT, J, I, MR, LT, JCSX, JX

      PARAMETER (MT = 127,MR = 90)
      COMMON /STAT1 / CSX(10), CS(90,127,10), CNIU(60,41,2), W(18400),
     :                JCSX
      COMMON /STAT3 / F(4,4,127,90,10),TFM(127),FRM(90),
     :                FKLM(6),FKHM(6)

*
* Calculate a bicubic spline interpolation fit for the temperature
* and density opacity fit.  Do not stop if the input lies outside the
* array but rather use the previous result.
*
* Check that we are in the table.
*
      IF((TF.LT.TFM(1)).OR.(TF.GE.TFM(MT)).OR.
     &   (FR.LT.FRM(1)).OR.(FR.GE.FRM(MR))) THEN
         FKL = FKLO
         FKH = FKHO
         WRITE(6,100) TF,FR,FKL,FKH
      ELSE
         FKLO = FKL
         FKHO = FKH
*
* Find interval in which target point lies.
*
         I = 1 + (MT-1)*(TF-TFM(1))/(TFM(MT)-TFM(1))
         J = 1 + (MR-1)*(FR-FRM(1))/(FRM(MR)-FRM(1))
         DT = TF-TFM(I)
         DR = FR-FRM(J)
*
* Evaluate the splines.
*
         FKL = F(1,1,I,J,JX) + DR*(F(1,2,I,J,JX)
     &   + DR*(F(1,3,I,J,JX) + DR*F(1,4,I,J,JX)))
     &   + DT*(F(2,1,I,J,JX) + DR*(F(2,2,I,J,JX)
     &   + DR*(F(2,3,I,J,JX) + DR*F(2,4,I,J,JX)))
     &   + DT*(F(3,1,I,J,JX) + DR*(F(3,2,I,J,JX)
     &   + DR*(F(3,3,I,J,JX) + DR*F(3,4,I,J,JX)))
     &   + DT*(F(4,1,I,J,JX) + DR*(F(4,2,I,J,JX)
     &   + DR*(F(4,3,I,J,JX) + DR*F(4,4,I,J,JX))))))
*
         FKH = F(1,1,I,J,JX+1) + DR*(F(1,2,I,J,JX+1)
     &   + DR*(F(1,3,I,J,JX+1) + DR*F(1,4,I,J,JX+1)))
     &   + DT*(F(2,1,I,J,JX+1) + DR*(F(2,2,I,J,JX+1)
     &   + DR*(F(2,3,I,J,JX+1) + DR*F(2,4,I,J,JX+1)))
     &   + DT*(F(3,1,I,J,JX+1) + DR*(F(3,2,I,J,JX+1)
     &   + DR*(F(3,3,I,J,JX+1) + DR*F(3,4,I,J,JX+1)))
     &   + DT*(F(4,1,I,J,JX+1) + DR*(F(4,2,I,J,JX+1)
     &   + DR*(F(4,3,I,J,JX+1) + DR*F(4,4,I,J,JX+1))))))
      ENDIF
  100 FORMAT ('OPACITY OUT OF RANGE',' TF FR FKL FKH ',4F9.4)
      RETURN
      END


