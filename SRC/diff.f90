!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   calculate the df(x)/dx using nine point Richardson’s extrapolation formula  !   
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PURE FUNCTION DIFF(P,DX)
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  
  REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: P
  REAL(KIND=DBL), INTENT(IN) :: DX
  REAL(KIND=DBL), DIMENSION(1:SIZE(P)) :: DIFF

  INTEGER :: I, N

  N = SIZE(P)
  DIFF(1) = (P(N-1)-8.0D0*P(N)+8.0D0*P(2)-P(3))/(12.0D0*DX)
  DIFF(2) = (P(N)-8.0D0*P(1)+8.0D0*P(3)-P(4))/(12.0D0*DX)
  ! 5 point formula to calculate first order derivative at edges
  FORALL(I=3:4)
    DIFF(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END FORALL
  FORALL(I=N-3:N-2)
    DIFF(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END FORALL
  ! 9 point formula to calculate first order derivative
  FORALL(I=5:N-4)
    DIFF(I) = (3.0D0*P(I-4)-32.0D0*P(I-3)+168.0D0*P(I-2)-672.0D0*P(I-1) +&
              672.0D0*P(I+1)-168.0D0*P(I+2)+32.0D0*P(I+3)-3.0D0*P(I+4))/&
              (840.0D0*DX)
  END FORALL
  DIFF(N-1) = (P(N-3)-8.0D0*P(N-2)+8.0D0*P(N)-P(1))/(12.0D0*DX)
  DIFF(N) = (P(N-2)-8.0D0*P(N-1)+8.0D0*P(1)-P(2))/(12.0D0*DX)
END FUNCTION DIFF
