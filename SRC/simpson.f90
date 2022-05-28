!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            calculate the integration by Simpson’s 1/3 rule             !   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PURE FUNCTION SIMPSON(F, DX)
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: F
  REAL(KIND=DBL), INTENT(IN) :: DX
  REAL(KIND=DBL) :: SIMPSON
  REAL(KIND=DBL) :: F1, F2
  INTEGER :: I, N
  N = SIZE(F)  !N is even
  F1 = F(2) + F(N-4)
  F2 = F(3)
  DO I = 4, N-6, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  SIMPSON = DX*(F(1) + 4.0D0*F1 + 2.0D0*F2 + F(N-3))/3.0D0+&
         3.0D0*DX*(F(N-3) + 3.0D0*(F(N-2)+F(N-1)) + F(N))/8.0D0
END FUNCTION SIMPSON
