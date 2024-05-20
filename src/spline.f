C     Calculate the coefficients of a 1-D cubic spline:
C
C     Forsythe, Malcolm, Moler, Computer Methods for Mathematical
C     Computations, Prentice-Hall, 1977, p.76
      SUBROUTINE SPLINE(N,X,F)

      IMPLICIT NONE

      REAL*8 F, X, T
      INTEGER I, N, IB

      DIMENSION X(*),F(4,*)
*
      F(2,N) = 0D0
      F(3,N) = 0D0
      F(4,N) = 0D0
*
* Set up a tridiagonal system for A*y=B where y(i) are the second
* derivatives at the knots.
* f(2,i) are the diagonal elements of A
* f(4,i) are the off-diagonal elements of A
* f(3,i) are the B elements/3, and will become c/3 upon solution
*
      F(4,1) = X(2)-X(1)
      F(3,2) = (F(1,2) - F(1,1))/F(4,1)
      DO I = 2,N-1
         F(4,I) = X(I+1) - X(I)
         F(2,I) = 2D0*(F(4,I-1) + F(4,I))
         F(3,I+1) = (F(1,I+1) - F(1,I))/F(4,I)
         F(3,I) = F(3,I+1) - F(3,I)
      ENDDO
*
* Boundaries.
*
      F(2,2) = F(4,1) + 2D0*F(4,2)
      F(3,2) = F(3,2)*F(4,2)/(F(4,1) + F(4,2))
*
      F(2,N-1) = 2D0*F(4,N-2) + F(4,N-1)
      F(3,N-1) = F(3,N-1)*F(4,N-2)/(F(4,N-1) + F(4,N-2))
*
* Forward elimination.
*
      T = F(4,2)/F(2,2)
      F(2,3) = F(2,3) - T*(F(4,2) - F(4,1))
      F(3,3) = F(3,3) - T*F(3,2)
      DO I = 4,N-2
         T = F(4,I-1)/F(2,I-1)
         F(2,I) = F(2,I)-T*F(4,I-1)
         F(3,I) = F(3,I)-T*F(3,I-1)
      ENDDO
      T = (F(4,N-1) - F(4,N-1))/F(2,N-2)
      F(2,N-1) = F(2,N-1) - T*F(4,N-2)
      F(3,N-1) = F(3,N-1) - T*F(3,N-2)
*
* Back substitution.
*
      F(3,N-1) = F(3,N-1)/F(2,N-1)
      DO IB = 1,N-4
         I = N-1-IB
         F(3,I) = (F(3,I) - F(4,I)*F(3,I+1))/F(2,I)
      ENDDO
      F(3,2) = (F(3,2) - (F(4,2) - F(4,1))*F(3,3))/F(2,2)
*
* Reset d array to step size.
*
      F(4,1) = X(2) - X(1)
      F(4,N-1) = X(N) - X(N-1)
*
* Set f(3,1) for not-a-knot.
*
      F(3,1) = (F(3,2)*(F(4,1) + F(4,2)) - F(3,3)*F(4,1))/F(4,2)
      F(3,N) = F(3,N-1) + (F(3,N-1) - F(3,N-2))*F(4,N-1)/F(4,N-2)
*
* Compute the polynomial coefficients.
*
      DO I = 1,N-1
         F(2,I) = (F(1,I+1) - F(1,I))/F(4,I) - F(4,I)*(F(3,I+1) + 2.0D0*F(3,I))
         F(4,I) = (F(3,I+1) - F(3,I))/F(4,I)
         F(3,I) = 3D0*F(3,I)
         F(4,I) = F(4,I)
      ENDDO
      RETURN
      END
