C Compute yields from surface output in BS
      IMPLICIT NONE
      REAL*8 yield(50), age, crap(5)
      REAL*8 compos(50), masslost
      INTEGER I,J,nmod
      OPEN(1,file="input",status="old")
      OPEN(2,file="output",status="new")
      OPEN(3,file="initialcompos",status="old")
 117  format (I6,1P,E16.9,54(1X,E13.6))
      DO I=1,50
         YIELD(I) = 0d0
      END DO
      READ(3,117) NMOD,AGE,(compos(J),J=1,50),(crap(i),i=1,3)
      READ(1,117) NMOD,AGE,(yield(J),J=1,50),(crap(i),i=1,2),masslost
      DO J=1,50
         YIELD(J) = YIELD(J) + COMPOS(J)*masslost
      END DO
      write(2,117) NMOD,AGE,(yield(j),j=1,50),(crap(i),i=1,2),masslost
      END
