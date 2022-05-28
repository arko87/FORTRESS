! File name: cgpe1d.f90 
! Title: FORTRESS: FORTRAN programs to solve coupled
! Gross-Pitaevskii equations for spin-orbit coupled spin-f
! Bose-Einstein condensate with spin f = 1 or 2
!
! Authors: Paramjeet Banger, Pardeep Kaur, Arko Roy, & Sandeep Gautam
!
! Fortran code to solve 1D coupled Gross-Pitaevskii equations (CGPEs)
! for spin-2(1) BEC using time-splitting spectral method.       
!
!----------------------------------------------------------------------------
! PI = 4.0*arctan(1), CI = \sqrt{-1}                                                          
! NITER = Number of total iterations 
! NSTP  = Number of iterations after which transeint results are written
! OPENMP_THREADS = Number of threads chosen for OpenMP parallelization
! FFTW_THREADS = Number of FFTW threads
! NX = Number of space grid points
! DX = Spatial step size
! DT = Temporal step size
! LX = NX*DX = 1D spatial domain
! AMU = atomic mass unit in Kg
! HBAR = reduced Planck's constat in SI units
! CDT = Complex varible which will be defined as '-idt' and 'dt' in
!       imaginary and realtime propagations respectively 
!-----------------------------------------------------------------------------

MODULE  DOUBLE_PRECISION
  INTEGER, PARAMETER :: DBL = KIND(0.0D0)
END MODULE  DOUBLE_PRECISION

MODULE BASIC_DATA
  USE DOUBLE_PRECISION
  REAL(KIND=DBL), PARAMETER :: PI = 4.0D0*DATAN(1.0D0)
  COMPLEX(KIND=DBL), PARAMETER :: CI = CMPLX(0.0D0,1.0D0,KIND=DBL)
  REAL(KIND=DBL), PARAMETER :: AMU = 1.66D-27, HBAR = 1.054560653D-34, AU = 0.529177208D-10
  REAL(KIND=DBL), PARAMETER :: EPS = 1.0D-40
  REAL(KIND=DBL), PARAMETER :: TOL_NR = 1.0D-6
  INTEGER :: NITER
  INTEGER :: OPENMP_THREADS, FFTW_THREADS
  INTEGER :: NX, NX2, NXX
  REAL(KIND=DBL) :: DX, DT
  REAL(KIND=DBL) :: LX 
  INTEGER :: STP, NSTP
  COMPLEX(KIND=DBL) :: CDT
END MODULE BASIC_DATA

MODULE CGPE_DATA
  USE BASIC_DATA, ONLY : NX 
  USE DOUBLE_PRECISION
!------------------------------------------------------------------------------
  !Scattering length values for Na-23-Antiferromagnetic
  !REAL(KIND=DBL), PARAMETER :: M = (23.0D0*AMU)
  !REAL(KIND=DBL), PARAMETER :: A0 = 34.9D0*AU, A2 = 45.8D0*AU, A4 = 64.5D0*AU
  !Scattering length values for Rb-83-Ferromagnetic
  !REAL(KIND=DBL), PARAMETER :: M = (83.0D0*AMU) 
  !REAL(KIND=DBL), PARAMETER :: A0 = 83.0D0*AU, A2 = 82.0D0*AU, A4 = 81.0D0*AU
  !Scattering length values for Rb-87-Cyclic
  !REAL(KIND=DBL), PARAMETER :: M = (87.0D0*AMU) 
  !REAL(KIND=DBL), PARAMETER :: A0 = 87.93D0*AU, A2 = 91.28D0*AU, A4 = (99.18D0)*AU
!------------------------------------------------------------------------------
  REAL(KIND=DBL) :: M
  REAL(KIND=DBL) :: A0, A2, A4
  REAL(KIND=DBL) :: NUX, NUY, NUZ,&
                    ALPHAX, ALPHAY, ALPHAZ
  INTEGER :: NATOMS, SPIN 
  INTEGER, PARAMETER :: NTRIAL = 1000
  REAL(KIND=DBL) :: TOL
  INTEGER :: SWITCH_IM, OPTION_FPC
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: X, X2
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: KX
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: V, R2
  COMPLEX(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: PHI, PHIF, PHI_OLD
  REAL(KIND=DBL), DIMENSION(0:2) :: TAU
  REAL(KIND=DBL) :: AOSC, OMEGAM, MAG, N1, N2, N3, N4, N5
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: SX
  COMPLEX(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: SY
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: SZ
END MODULE CGPE_DATA

MODULE FFTW3
  USE, INTRINSIC :: ISO_C_BINDING
  INCLUDE 'fftw3.f03'
END MODULE FFTW3
 
MODULE SOC_DATA
  USE DOUBLE_PRECISION
  REAL(KIND=DBL) :: GAMMAX
  INTEGER :: SWITCH_SOC 
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: A
END MODULE SOC_DATA

MODULE FFTW_DATA
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : NX
  USE FFTW3, ONLY : C_DOUBLE_COMPLEX, C_PTR, C_INT
  COMPLEX(KIND=C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:) :: FFTFX, FFTBX
  TYPE(C_PTR) :: PLANFX, PLANBX
  INTEGER(KIND=C_INT) :: THREADS_INIT
END MODULE FFTW_DATA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                      Main Program - CGPE1D                                       ! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PROGRAM CGPE1D
  USE BASIC_DATA
  USE CGPE_DATA
  USE SOC_DATA
  USE FFTW3, ONLY : C_PTR, FFTW_INIT_THREADS, FFTW_PLAN_WITH_NTHREADS, FFTW_CLEANUP_THREADS,&
                    FFTW_PLAN_DFT_1D, FFTW_EXECUTE_DFT, FFTW_DESTROY_PLAN
  USE OMP_LIB
  USE FFTW_DATA
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE ALLOCATE_MEM()
      IMPLICIT NONE
    END SUBROUTINE ALLOCATE_MEM

    SUBROUTINE DEALLOCATE_MEM()
      IMPLICIT NONE
    END SUBROUTINE DEALLOCATE_MEM
 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
    END SUBROUTINE INITIALIZE
    
    SUBROUTINE KE()
      IMPLICIT NONE
    END SUBROUTINE KE
   
    SUBROUTINE SOC()
     IMPLICIT NONE
    END SUBROUTINE SOC

    SUBROUTINE SP()
      IMPLICIT NONE
    END SUBROUTINE SP
  
    SUBROUTINE SE()
      IMPLICIT NONE
    END SUBROUTINE SE
 
    SUBROUTINE NORMC(PHI, NORM)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: NORM
    END SUBROUTINE NORMC
    
    SUBROUTINE NORMT(PHI, TNORM)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(INOUT) :: PHI
      REAL(KIND=DBL), INTENT(OUT) :: TNORM
    END SUBROUTINE NORMT 

    SUBROUTINE RAD(PHI, RMS)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:, :), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: RMS
    END SUBROUTINE RAD

    SUBROUTINE ENERGY(PHI, MU, EN, MZ)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:, :), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: MU
      REAL(KIND=DBL), INTENT(OUT) :: EN, MZ
    END SUBROUTINE ENERGY

    SUBROUTINE NEWTON(X, N1, N2, N3, N4, N5, MAG, ERR_MSG)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), INTENT(IN) :: N1, N2, N3, N4, N5, MAG
      INTEGER, INTENT(INOUT) :: ERR_MSG
      REAL(KIND=DBL), DIMENSION(1:2), INTENT(INOUT) :: X
    END SUBROUTINE NEWTON

    SUBROUTINE FFT()
      IMPLICIT NONE
    END SUBROUTINE FFT

    SUBROUTINE BFT()
      IMPLICIT NONE
    END SUBROUTINE BFT
    
    SUBROUTINE CREATE_PLANS()
      IMPLICIT NONE
    END SUBROUTINE CREATE_PLANS

    SUBROUTINE DESTROY_PLANS()
      IMPLICIT NONE
    END SUBROUTINE DESTROY_PLANS
  END INTERFACE
  
  INTEGER :: I, K, MS, IOSTAT, II, ERR_MSG
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: NORM, RMS, MU, SIGMA
  REAL(KIND=DBL) :: EN, EN_OLD, MZ, START, FINISH
  REAL(KIND=DBL) :: TNORM
  REAl(KIND=DBL), DIMENSION(1:2) :: SIG
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: TMP
  REAL(KIND=DBL) :: CONV = 0.0D0, CONV_EN = 0.0D0
  CHARACTER(LEN=10) :: DATE, TIME, ZONE
  INTEGER, DIMENSION(1:8) :: VALUES
  CHARACTER(LEN=1) :: SPIN_ID
  LOGICAL :: FEXIST

  !$ START = OMP_GET_WTIME()
  CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)

  OPEN (100, FILE='err_msg.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
              ACTION='WRITE')

  INQUIRE (FILE='./INPUT/input1D.dat', EXIST = FEXIST)
  IF(FEXIST)THEN
     OPEN (10, FILE='./INPUT/input1D.dat', STATUS='OLD',FORM='FORMATTED',&
              ACTION='READ')
  ELSE
     WRITE(100,*) "'input1D.dat' does not exist."
     PRINT*, "'input1D.dat' does not exist."
     STOP
  END IF

  II = 0
  ERR_MSG = 0
  SWITCH_SOC = 0
  OPTION_FPC = 4

  READ(10,*) SPIN
  READ(10,*) NITER
  READ(10,*) STP, NSTP
  READ(10,*) OPENMP_THREADS
  READ(10,*) NX
  READ(10,*) DX, DT
  READ(10,*) M
  IF(SPIN.EQ.1)THEN
     READ(10,*) A0, A2
  ELSE
     READ(10,*) A0, A2, A4
  END IF
  READ(10,*) NATOMS
  READ(10,*) NUX, NUY, NUZ
  READ(10,*) TOL
  READ(10,*) SWITCH_IM
  READ(10,*) GAMMAX
  IF(ABS(GAMMAX).GT.EPS) SWITCH_SOC = 1  
  IF(SWITCH_SOC.EQ.0)THEN
     READ(10,*) OPTION_FPC
  ELSE
     READ(10,*)
  END IF
  IF(SWITCH_IM.EQ.1.AND.SWITCH_SOC.EQ.0)THEN
     READ(10,*) MAG
  ELSE
     READ(10,*)
  END IF
  CLOSE(10) 

  FFTW_THREADS = OPENMP_THREADS
  NX2 = NX/2+1
  NXX = NX-1

  ALPHAX = NUX/NUX
  ALPHAY = NUY/NUX
  ALPHAZ = NUZ/NUX

  LX = NX*DX

  ALLOCATE(NORM(1:2*SPIN+1), RMS(1:2*SPIN+1), MU(1:2*SPIN+1), SIGMA(1:2*SPIN+1))
  ALLOCATE(TMP(1:4*SPIN+2))
  CALL ALLOCATE_MEM()

  !$ CALL OMP_set_num_threads(OPENMP_THREADS)
  !$OMP PARALLEL
    !$OMP SINGLE
    !$ WRITE(*,*) 'MAXIMUM_THREAD_NUMBER_SET = ', OMP_get_num_threads()
    !$OMP END SINGLE
  !$OMP END PARALLEL

  IF(SPIN.EQ.1)THEN 
     WRITE(SPIN_ID, '(I1)') 1
  ELSE 
     WRITE(SPIN_ID, '(I1)') 2
  END IF

  SELECT CASE(SWITCH_IM)
    CASE(1)
      CDT = -CI*DT
      PHI_OLD = CMPLX(0.0D0,0.0D0,KIND=DBL)
      EN_OLD = 0.0D0
      OPEN (2 ,FILE='./OUTPUT/file1_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
      OPEN (3, FILE='./OUTPUT/file2_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
      OPEN (4, FILE='./OUTPUT/file3_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
      OPEN (8, FILE='./OUTPUT/solution_file_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
      OPEN (9, FILE='./OUTPUT/convergence.dat', STATUS='UNKNOWN', FORM='FORMATTED',&
               ACTION='WRITE')
    CASE(0)
      CDT = CMPLX(DT,0.0D0,KIND=DBL)
      OPEN (1, FILE='./INPUT/initial_sol_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               IOSTAT=iostat, STATUS='OLD',FORM='FORMATTED', ACTION='READ')
      OPEN (2, FILE='./OUTPUT/file1_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
      OPEN (3, FILE='./OUTPUT/file2_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
      OPEN (4, FILE='./OUTPUT/file3_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
      OPEN (8, FILE='./OUTPUT/solution_file_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
               STATUS='UNKNOWN',FORM='FORMATTED', ACTION='WRITE')
  END SELECT

  WRITE(2,1007) VALUES(3),VALUES(2),VALUES(1),VALUES(5),VALUES(6),VALUES(7)

  OMEGAM = 2.0D0*PI*NUX
  AOSC = SQRT(HBAR/(M*AMU*OMEGAM))

  SELECT CASE(SPIN)
    CASE(1)
      TAU(0) = 2.0D0*(SQRT(ALPHAY*ALPHAZ)/(2.0D0*PI))*(4.0D0*PI*REAL(NATOMS)*(A0*AU+2.0D0*A2*AU)/(3.0D0*AOSC))
      TAU(1) = 2.0D0*(SQRT(ALPHAY*ALPHAZ)/(2.0D0*PI))*(4.0D0*PI*(REAL(NATOMS))*(A2*AU-A0*AU)/(3.0D0*AOSC))
    CASE(2)
      TAU(0) = 2.0D0*((SQRT(ALPHAY*ALPHAZ)/(2.0D0*PI))*4.0D0*PI*DBLE(NATOMS)*(4.0D0*A2*AU+3.0D0*A4*AU)/(7.0D0*AOSC))
      TAU(1) = 2.0D0*((SQRT(ALPHAY*ALPHAZ)/(2.0D0*PI))*4.0D0*PI*DBLE(NATOMS)*(A4*AU-A2*AU)/(7.0D0*AOSC))
      TAU(2) = 2.0D0*((SQRT(ALPHAY*ALPHAZ)/(2.0D0*PI))*4.0D0*PI*DBLE(NATOMS)*(7.0D0*A0*AU-10.0D0*A2*AU+3.0D0*A4*AU)/(7.0D0*AOSC))
  END SELECT

  WRITE(2,*)
  WRITE(2,899) OPENMP_THREADS, FFTW_THREADS
  WRITE(2,900) SWITCH_IM, OPTION_FPC, SWITCH_SOC, GAMMAX, TOL
  WRITE(2,901) ALPHAX, ALPHAY,ALPHAZ
  WRITE(2,*)
  WRITE(2,902) NX
  WRITE(2,903) NITER, NSTP, STP
  WRITE(2,905) DX
  WRITE(2,906) DT, MAG
  WRITE(2,*) 
  WRITE(2,907) AOSC
  IF (SPIN .EQ. 1) THEN
     WRITE(2,904) A0, A2
     WRITE(2,908) TAU(0)/2.0D0, TAU(1)/2.0D0
  ELSE 
     WRITE(2,911) A0, A2, A4
     WRITE(2,912) TAU(0)/2.0D0, TAU(1)/2.0D0, TAU(2)/2.0D0 
  END IF
  WRITE(2,*)

  899 FORMAT(' OPENMP_THREADS = ',I3,', FFTW_THREADS = ', I3)
  900 FORMAT(' SWITCH_IM = ',I2,', OPTION_FPC = ', I2, ', SWITCH_SOC = ', I2, ', GAMMAX = ', F6.2, ', TOL = ', ES9.2)
  901 FORMAT(' Anisotropy ALPHAX = ',F12.6,', ALPHAY = ',F12.6,', ALPHAZ = ',F12.6)
  902 FORMAT(' No. of space steps NX = ',I8)
  903 FORMAT(' No. of time steps : NITER = ',I9,', NSTP = ',I9,', STP = ',I7)
  904 FORMAT(' A0 = ', F12.6, ', A2 = ', F12.6, ',  (all in units of Bohr radius) ')
  905 FORMAT(' Space Step DX = ',F10.6)
  906 FORMAT(' Time Step  DT = ', F10.6, ', MAG =',F12.6)
  907 FORMAT(' AOSC = ', E12.4 )
  908 FORMAT(' TAU(0) = ', F12.6,', TAU(1) = ', F12.6,' (in dimensionless units)')
  911 FORMAT(' A0 = ', F12.6, ', A2 = ', F12.6, ', A4 = ', F12.6, ' (all in units of Bohr radius) ')
  912 FORMAT(' TAU(0) = ', F12.6,', TAU(1) = ', F12.6,', TAU(2) = ', F12.6, ' (in dimensionless units)')

  CALL INITIALIZE()

  IF((SWITCH_IM.NE.1))THEN
    IF(IOSTAT.NE.0) THEN
      PRINT*, 'Check if the "initial_sol_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" is available?'
      PRINT*,'if not, user can read the "solution_file_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" file as an initial solution.'
      WRITE(100,*)'Check if the "initial_sol_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" is available?'
      WRITE(100,*)'if not, user can read the "solution_file_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" file as an initial solution.'
      STOP
    END IF
    DO I = 1, NX
       READ(1,*) X(I), (TMP(K), K = 1, 4*SPIN+2)
       SELECT CASE(SPIN)
         CASE(1)
           PHI(I,1) = SQRT(TMP(1))*EXP(CI*TMP(2))
           PHI(I,2) = SQRT(TMP(3))*EXP(CI*TMP(4))
           PHI(I,3) = SQRT(TMP(5))*EXP(CI*TMP(6))
         CASE(2)
           PHI(I,1) = SQRT(TMP(1))*EXP(CI*TMP(2))
           PHI(I,2) = SQRT(TMP(3))*EXP(CI*TMP(4))
           PHI(I,3) = SQRT(TMP(5))*EXP(CI*TMP(6))
           PHI(I,4) = SQRT(TMP(7))*EXP(CI*TMP(8))
           PHI(I,5) = SQRT(TMP(9))*EXP(CI*TMP(10))
       END SELECT
    END DO
  END IF
  
  IF (SWITCH_IM .EQ. 1) CALL NORMT(PHI,TNORM)

  MS = 2*SPIN+1

  CALL NORMC(PHI, NORM)
  CALL RAD(PHI, RMS)
  CALL ENERGY(PHI, MU, EN, MZ)

  SELECT CASE(SPIN)
    CASE(1)
      WRITE (2, 1001)
      WRITE (2, 1002)
      WRITE (2, 1001)
      WRITE (2, 1003) SUM(NORM), EN/2.0D0, (MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,K)), K = 1, MS)
    CASE(2)
      WRITE (2, 1011)
      WRITE (2, 1012)
      WRITE (2, 1011)
      WRITE (2, 1013) SUM(NORM), EN/2.0D0, (MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,K)), K = 1, MS)
  END SELECT

  1001 FORMAT (23X, 76('-'))  
  1002 FORMAT (23X, 'Norm', 7X,'EN', 7X, 'MU1', 7X, 'MU2', 7X, 'MU3', 6X,&
                    'phi1', 6X, 'phi2', 6X,'phi3')
  1003 FORMAT ('Initial: ', 8X, 1F11.4, 4F10.5, 3F10.5)
  1011 FORMAT (22X, 117('-'))
  1012 FORMAT (23X, 'Norm', 7X,'EN', 7X, 'MU1', 7X, 'MU2', 7X, 'MU3', 7X, 'MU4', 7X, 'MU5', 7X,&
                    'phi1', 6X, 'phi2', 6X,'phi3',5X, 'phi4', 6X,'phi5')
  1013 FORMAT ('Initial : ', 7X, 1F11.4, 6F10.5, 5F10.5)

  !Initializing SIG(1) and SIG(2) 
  SIG(1) = 20.0D0
  SIG(2) = 20.0D0

  CALL CREATE_PLANS()

  WRITE(3,'(F14.6,6F18.10)') 2.0D0*II*DT, EN/2.0D0, (RMS(K), K = 1, MS)
  WRITE(4,'(F14.6,7F18.10)') 2.0D0*II*DT, (NORM(K), K = 1, MS), SUM(NORM), MZ

  DO II = 1, NITER
     CALL FFT()
     CALL KE()
     CALL SOC()
     CALL BFT()
     CALL SE()
     CALL SP()
     
     CALL NORMC(PHI, NORM)
     
     IF(SWITCH_IM.EQ.1)THEN
       IF(SWITCH_SOC.EQ.0)THEN
         SELECT CASE (SPIN) 
           CASE(1)
             SIGMA(2) = SQRT(1.0D0 - MAG**2)/(SQRT(NORM(2)+SQRT(4.0D0*&
                      (1.0D0-MAG**2)*NORM(1)*NORM(3)+MAG**2*NORM(2)**2))+EPS)
             SIGMA(1) = SQRT(1.0D0+MAG-SIGMA(2)**2*NORM(2))/(SQRT(2.0D0*NORM(1))+EPS)
             SIGMA(3) = SQRT(1.0D0-MAG-SIGMA(2)**2*NORM(2))/(SQRT(2.0D0*NORM(3))+EPS)
           CASE(2)
             CALL NEWTON(SIG,NORM(1),NORM(2),NORM(3),NORM(4),NORM(5),MAG,ERR_MSG)
             IF(ERR_MSG.NE.0) STOP
             !Projection operators
             SIGMA(1) = EXP((SIG(1) + 2.0D0*SIG(2))*DT)
             SIGMA(2) = EXP((SIG(1) + SIG(2))*DT)
             SIGMA(3) = EXP(SIG(1)*DT)
             SIGMA(4) = EXP((SIG(1) - SIG(2))*DT)
             SIGMA(5) = EXP((SIG(1) - 2.0D0*SIG(2))*DT) 
         END SELECT         

         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,K)
         DO K = 1, MS
            DO I = 1, NX
               PHI(I,K) = SIGMA(K)*PHI(I,K)
            END DO
         END DO
         !$OMP END PARALLEL DO        

         CALL NORMT(PHI, TNORM)
       ELSE
         CALL NORMT(PHI, TNORM)
       END IF
     END IF
 
     IF(MOD(II,STP).EQ.0)THEN 
       CALL ENERGY(PHI, MU, EN, MZ) 
       CALL RAD(PHI, RMS)
    
       IF(SWITCH_IM.EQ.1)THEN
          !$OMP PARALLEL WORKSHARE
          CONV = MAXVAL(ABS(PHI(:,:)-PHI_OLD(:,:)))/(STP*2.0D0*DT)
          CONV_EN = ABS(EN-EN_OLD)/(STP*2.0D0*DT)
          !$OMP END PARALLEL WORKSHARE
          WRITE(9,*) 2.0D0*II*DT, CONV, CONV_EN
          IF(CONV.LE.TOL.OR.CONV_EN.LE.TOL) EXIT
       END IF
       !-------------------------------------------------------------------------------------------!
       ! In file2_*_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat, time is written in first              !
       ! column energy in second column the rms sizes in the last (2*SPIN+1) columns.              !            
       !                                                                                           !
       ! In file3_*_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat, time is written in first              !
       ! column, norms of the individual components in the next (2*SPIN+1) columns and             !               
       ! magnetization in last column.                                                             !
       !                                                                                           ! 
       ! In tmp_solution_file_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat, first column has x             !
       ! variable, density of first component followed by its phase in second and third            !
       ! columns, respectively, then the same for all components respectively.                     !
       !                                                                                           !
       ! In convergence.dat, time, wavefunction convergence and energy convergence                 ! 
       ! parameter are written in  first and second columns, respectively.                         ! 
       !-------------------------------------------------------------------------------------------!

       WRITE(3,'(F14.6,6F18.10)') 2.0D0*II*DT, EN/2.0D0, (RMS(K), K=1,MS)
       WRITE(4,'(F14.6,7F18.10)') 2.0D0*II*DT, (NORM(K), K=1,MS), SUM(NORM), MZ
     END IF    
     EN_OLD = EN
     PHI_OLD(:,:) = PHI(:,:)

     IF(MOD(II,NSTP).EQ.0)THEN
        WRITE(2, 1005) SUM(NORM), EN/2.0D0, (MU(K)/2.0D0,K=1,MS), (ABS(PHI(NX2,K)), K=1,MS)
        OPEN(7, FILE='./OUTPUT/tmp_solution_file_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
                STATUS='replace', FORM='FORMATTED', ACTION='WRITE')
        DO I = 1, NX
           WRITE(7,'(F11.5,10F14.6)') X(I), (ABS(PHI(I,K))**2, ATAN2(AIMAG(PHI(I,K)), REAL(PHI(I,K))), K=1,MS)
        END DO
        CLOSE(7)
     END IF
  END DO   

  CALL DESTROY_PLANS()
  CALL ENERGY(PHI, MU, EN, MZ)
  CALL RAD(PHI, RMS)

  IF (SPIN.EQ.1) THEN
      WRITE (2, 1001)
  ELSE
     WRITE (2, 1011)
  END IF

  WRITE (2, 1006) II-1, SUM(NORM), EN/2.0D0, (MU(K)/2.0D0, K = 1,MS), (ABS(PHI(NX2,K)), K = 1,MS)
  1005 FORMAT('After NSTP iter.    :', 1F7.4, 6F10.5, 5F10.5)
  1006 FORMAT('After ', I8 ,' iter.:', 1F7.4, 6F10.5, 5F10.5)
 
  IF (SPIN.EQ.1) THEN
     WRITE (2, 1001)
  ELSE
     WRITE (2, 1011)
  END IF
  WRITE(2,*)
 
  !------------------------------------------------------------------------------------------------
  ! In solution_file.dat, first column has x variable, density of first component followed by     !
  ! its phase in second and third columns,  respectively, then the same for all components.       !     
  !------------------------------------------------------------------------------------------------

  DO I = 1, NX
     WRITE(8,'(F11.5,10F18.10)') X(I), (ABS(PHI(I,K))**2, ATAN2(AIMAG(PHI(I,K)),REAL(PHI(I,K))), K = 1,MS)
  END DO
  CLOSE(8)
  CALL DEALLOCATE_MEM()
  PRINT *,'Program execution complete'

  !$ FINISH = OMP_GET_WTIME()
  CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
  WRITE(2,1008) VALUES(3),VALUES(2),VALUES(1),VALUES(5),VALUES(6),VALUES(7)
  !$ WRITE(2,1009) (FINISH-START)
  1007 FORMAT("Started on DAY/MONTH/YEAR, HOURS:MINUTES:SECONDS = ",I2,'/',I2,'/',I4,',',I3,':',I2,':',I2)
  1008 FORMAT("Ended on DAY/MONTH/YEAR, HOURS:MINUTES:SECONDS = ",I2,'/',I2,'/',I4,',',I3,':',I2,':',I2)
  1009 FORMAT("Time elapsed = ",F18.6, ' in seconds')

  CLOSE(100)
END PROGRAM CGPE1D

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Dynamic Memory Allocation to the variables and arrays                          !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE ALLOCATE_MEM()
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : NX
  USE CGPE_DATA, ONLY : SPIN, X, X2, V, R2, KX, PHI, PHIF, PHI_OLD, SX, SY, SZ
  USE SOC_DATA, ONLY : A
  USE FFTW_DATA, ONLY : FFTFX, FFTBX

  IMPLICIT NONE
  INTEGER :: MS, STAT

  MS = 2*SPIN+1
  ALLOCATE(X(1:NX),X2(1:NX),V(1:NX), R2(1:NX),KX(1:NX),STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(PHI(1:NX,1:MS),PHIF(1:NX,1:MS), PHI_OLD(1:NX,1:MS),STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(A(1:MS,1:MS), STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(SX(1:MS,1:MS),SY(1:MS,1:MS),SZ(1:MS), STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(FFTFX(NX), FFTBX(NX), STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
END SUBROUTINE ALLOCATE_MEM

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                       Memory Deallocation of the variables and arrays                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DEALLOCATE_MEM()
  USE DOUBLE_PRECISION
  USE CGPE_DATA, ONLY : X, X2, V, R2, KX, PHI, PHIF, PHI_OLD, SX, SY, SZ
  USE SOC_DATA, ONLY : A
  USE FFTW_DATA, ONLY : FFTFX, FFTBX

  IMPLICIT NONE
  INTEGER :: STAT

  DEALLOCATE(X, X2, V, R2, KX, STAT = STAT)
  DEALLOCATE(PHI, PHIF, PHI_OLD,STAT = STAT)
  DEALLOCATE(A)
  DEALLOCATE(SX,SY,SZ)
  DEALLOCATE(FFTFX, FFTBX)
END SUBROUTINE DEALLOCATE_MEM

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Defines the A, SX, SY, SZ, spatial and Fourier grids, trapping potential,           !
!              and initializes the component wave functions                                        !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INITIALIZE()
  USE BASIC_DATA, ONLY : CI, PI, NX, LX
  USE CGPE_DATA
  USE DOUBLE_PRECISION
  USE SOC_DATA, ONLY : A
  !USE OMP_LIB
                     
  IMPLICIT NONE
  REAL(KIND=DBL), DIMENSION(1:NX):: TMP1D, TFDENSITY
  INTEGER :: I 
  
  SELECT CASE(SPIN)
    CASE(1)
      A(1,1) = 0.5D0
      A(2,1) = -SQRT(0.5D0)
      A(3,1) = A(1,1)

      A(1,2) = A(2,1)
      A(2,2) = 0.0D0
      A(3,2) = -A(2,1)

      A(1,3) = A(1,1)
      A(2,3) = -A(2,1)
      A(3,3) = A(1,1)

      SX(:,:) = (SQRT(0.5D0))*RESHAPE([0.0D0,1.0D0,0.0D0,&
                1.0D0,0.0D0,1.0D0,&
                0.0D0,1.0D0,0.0D0],[3,3])
      SY(:,:) = (CI/SQRT(2.0D0))*RESHAPE([0.0D0,-1.0D0,0.0D0,&
                1.0D0,0.0D0,-1.0D0,&
                0.0D0,1.0D0,0.0D0],[3,3])
      SZ(:) = [1,0,-1]
    CASE(2)
      A(1,1) = 0.25D0
      A(2,1) = -0.5D0
      A(3,1) = SQRT(0.375D0)
      A(4,1) = A(2,1)
      A(5,1) = A(1,1)
 
      A(1,2) = A(2,1)
      A(2,2) = -A(2,1)
      A(3,2) = 0.0D0
      A(4,2) = A(2,1)
      A(5,2) = -A(2,1)
   
      A(1,3) = A(3,1)
      A(2,3) = A(3,2)
      A(3,3) = A(2,1)
      A(4,3) = A(3,2)
      A(5,3) = A(3,1)
   
      A(1,4) = A(2,1)
      A(2,4) = A(2,1)
      A(3,4) = A(3,2)
      A(4,4) = -A(2,1)
      A(5,4) = -A(2,1)

      A(1,5) = A(1,1)
      A(2,5) = -A(2,1)
      A(3,5) = A(3,1)
      A(4,5) = -A(2,1)
      A(5,5) = A(1,1)

      SX(:,:) = RESHAPE([0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,&
                1.0D0,0.0D0,SQRT(1.5D0),0.0D0,0.0D0,&
                0.0D0,SQRT(1.5D0),0.0D0,SQRT(1.5D0),0.0D0,&
                0.0D0,0.0D0,SQRT(1.5D0),0.0D0,1.0D0,&
                0.0D0,0.0D0,0.0D0,1.0D0,0.0D0],[5,5])
      SY(:,:) = CI*RESHAPE([0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,&
                   -1.0D0,0.0D0,SQRT(1.5D0),0.0D0,0.0D0,&
                   0.0D0,-SQRT(1.5D0),0.0D0,SQRT(1.5D0),0.0D0,&
                   0.0D0,0.0D0,-SQRT(1.5D0),0.0D0,1.0D0,&
                   0.0D0,0.0D0,0.0D0,-1.0D0,0.0D0],[5,5])
      SZ(:) = [2,1,0,-1,-2]
  END SELECT
  !$OMP PARALLEL SECTIONS   

  !$OMP SECTION
  DO I = 1, 1 + NX/2
     KX(I) = (I-1.0D0)*2.0D0*PI/LX
  END DO
 
  DO I = 1, NX/2 -1
     KX(I+1+NX/2) = -KX(1-I+NX/2)
  END DO

  !$OMP SECTION
  DO I = 1, NX
     X(I) = (-1.0D0 +2.0D0*(I-1.0D0)/NX)*LX/2.0D0
     X2(I) = X(I)*X(I)
  END DO
  !$OMP END PARALLEL SECTIONS
 
  !$OMP PARALLEL DO PRIVATE(I)
  DO I = 1, NX
     R2(I) = X2(I) 
     V(I) = ALPHAX*ALPHAX*X2(I)
     TMP1D(I) = ALPHAX*X2(I)
     SELECT CASE (SPIN)

       CASE(1)
         SELECT CASE (OPTION_FPC)
           CASE(1)!Initial guess wavefunctions for Ferromagnetic interactions
             IF(((3.0D0*(TAU(0)+TAU(2))/(4.0D0))**(0.67D0)- ALPHAX*ALPHAX*X2(I)).GE.0.0D0)THEN
               TFDENSITY(I) =  ((3.0D0*(TAU(0)+TAU(2))/(4.0D0))**(0.67D0)- ALPHAX*ALPHAX*X2(I))/((TAU(0)+TAU(2)))
             ELSE
               TFDENSITY(I) = 0.0D0
             END IF

             IF(((3.0D0*(TAU(0)+TAU(2))/(4.0D0))**(0.33D0)).GE.10.0D0)THEN
               PHI(I,1) = ((1.0D0+MAG)/2.0D0)*SQRT(TFDENSITY(I))
               PHI(I,2) = SQRT((1.0D0-MAG**2)/2.0D0)*SQRT(TFDENSITY(I))
               PHI(I,3) = ((1.0D0-MAG)/2.0D0)*SQRT(TFDENSITY(I))
             ELSE
               PHI(I,1) = ((1.0D0+MAG)/2.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
               PHI(I,2) = SQRT((1.0D0-MAG**2)/2.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
               PHI(I,3) = ((1.0D0-MAG)/2.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             END IF

           CASE(2)!Intitial guess wavefunctions for antiferromagnetic interactions 
             PHI(I,1) = (SQRT((1.0D0+MAG)/2.0D0))*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,2) = 0.0D0
             PHI(I,3) = (SQRT((1.0D0-MAG)/2.0D0))*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))

           CASE DEFAULT!For choosing Gaussian intial guess wavefunction
             PHI(I,1) = EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,2) = EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,3) = EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
         END SELECT
       CASE(2)
         SELECT CASE (OPTION_FPC)    
           CASE(1)!Initial guess wavefunctions for Ferromagnetic interactions
             IF(((3.0D0*(TAU(0)+4.0D0*TAU(1))/(4.0D0))**(0.67D0) - ALPHAX*ALPHAX*X2(I)).GE.0.0D0)THEN
               TFDENSITY(I) = ((3.0D0*(TAU(0)+4.0D0*TAU(1))/(4.0D0))**(0.67D0)-&
                                ALPHAX*ALPHAX*X2(I))/(TAU(0)+4.0D0*TAU(1))
             ELSE
               TFDENSITY(I) = 0.0D0
             END IF

             IF((3.0D0*(TAU(0)+4.0D0*TAU(1)/(4.0D0))**(0.33D0)).GE.10.0D0)THEN
               PHI(I,1) = (((2.0D0+MAG)**2)/16.0D0)*SQRT(TFDENSITY(I))
               PHI(I,2) = ((SQRT(4.0D0-MAG**2)*(2.0D0+MAG))/8.0D0)*SQRT(TFDENSITY(I))
               PHI(I,3) = ((SQRT(1.5D0)*(4.0D0-MAG**2))/8.0D0)*SQRT(TFDENSITY(I))
               PHI(I,4) = ((SQRT(4.0D0-MAG**2)*(2.0D0-MAG))/8.0D0)*SQRT(TFDENSITY(I))
               PHI(I,5) = (((2.0D0-MAG)**2)/16.0D0)*SQRT(TFDENSITY(I))
             ELSE
               PHI(I,1) = (((2.0D0+MAG)**2)/16.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
               PHI(I,2) = ((SQRT(4.0D0-MAG**2)*(2.0D0+MAG))/8.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
               PHI(I,3) = ((SQRT(1.5D0)*(4.0D0-MAG**2))/8.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
               PHI(I,4) = ((SQRT(4.0D0-MAG**2)*(2.0D0-MAG))/8.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
               PHI(I,5) = (((2.0D0-MAG)**2)/16.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             END IF
  
           CASE(2)!Intitial guess wavefunctions for polar interactions 
             PHI(I,1) = (SQRT(2.0D0+MAG)/2.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,2) = 0.0D0
             PHI(I,3) = 0.0D0
             PHI(I,4) = 0.0D0
             PHI(I,5) = (SQRT(2.0D0-MAG)/2.0D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))

           CASE(3)!Intitial guess wavefunctions for cyclic interactions 
             PHI(I,1) = (SQRT((1.0D0+MAG)/3.0D0))*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,2) = 0.0D0
             PHI(I,3) = 0.0D0
             PHI(I,4) = (SQRT((2.0D0-MAG)/3.0D0))*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,5) = 0.0D0

           CASE DEFAULT !For choosing Gaussian intial guess wavefunction
             PHI(I,1) =SQRT(0.2D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,2) =SQRT(0.2D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,3) =SQRT(0.2D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,4) =SQRT(0.2D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
             PHI(I,5) =SQRT(0.2D0)*EXP(-TMP1D(I)/2.0D0)/SQRT(SQRT(PI))
         END SELECT
     END SELECT
  END DO
 !$OMP END PARALLEL DO 
END SUBROUTINE INITIALIZE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Creates FFTW plans for forward and backward transforms                         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CREATE_PLANS()
  USE BASIC_DATA
  USE FFTW3
  USE FFTW_DATA

  IMPLICIT NONE    

  THREADS_INIT = FFTW_INIT_THREADS()
  CALL FFTW_PLAN_WITH_NTHREADS(FFTW_THREADS)
  PLANFX = FFTW_PLAN_DFT_1D(NX,FFTFX,FFTBX,FFTW_FORWARD,FFTW_ESTIMATE)
  PLANBX = FFTW_PLAN_DFT_1D(NX,FFTBX,FFTFX,FFTW_BACKWARD,FFTW_ESTIMATE)
END SUBROUTINE CREATE_PLANS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                    Destroys FFTW plans for forward and backward transforms                       !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DESTROY_PLANS()
  USE FFTW3
  USE FFTW_DATA

  IMPLICIT NONE

  CALL FFTW_DESTROY_PLAN(PLANFX)
  CALL FFTW_DESTROY_PLAN(PLANBX)
  CALL FFTW_CLEANUP_THREADS()
END SUBROUTINE DESTROY_PLANS
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      Solves the Hamiltonian corresponding to Kinetic Energy terms in Fourier space               !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KE()
  USE BASIC_DATA, ONLY : CI, NX, CDT
  USE CGPE_DATA, ONLY : KX, PHIF, SPIN
  !USE OMP_LIB

  IMPLICIT NONE
  INTEGER :: I, K, MS

  MS = 2*SPIN+1
  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,K)
  DO K = 1, MS
     DO I = 1, NX
        PHIF(I,K) = EXP(-CI*CDT*KX(I)*KX(I))*PHIF(I,K)
     END DO
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE KE
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    Solves the Hamiltonian corresponding to trapping potential and diagonal interaction terms     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SP()
  USE BASIC_DATA, ONLY : CI, NX, CDT
  USE CGPE_DATA, ONLY : V, TAU, PHI, SPIN
  USE DOUBLE_PRECISION
  !USE OMP_LIB
  
  IMPLICIT NONE

  REAL(KIND=DBL), DIMENSION(1:2*SPIN+1) :: TMP
  REAL(KIND=DBL), DIMENSION(1:NX) :: RHO, FZ
  INTEGER :: I, J, MS

  MS = 2*SPIN+1

  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO PRIVATE(I,TMP,J)
      DO I = 1, NX
         RHO(I) = SUM(ABS(PHI(I,1:MS))**2)
         TMP(1) = V(I) + TAU(0)*RHO(I) + TAU(1)*(ABS(PHI(I,1))**2 + ABS(PHI(I,2))**2&
                  -ABS(PHI(I,3))**2)
         TMP(2) = V(I) + TAU(0)*RHO(I) + TAU(1)*(ABS(PHI(I,1))**2 + ABS(PHI(I,3))**2)
         TMP(3) = V(I) + TAU(0)*RHO(I) + TAU(1)*(ABS(PHI(I,3))**2 + ABS(PHI(I,2))**2&
                  -ABS(PHI(I,1))**2)
         DO J = 1, MS
            PHI(I,J) = PHI(I,J)*EXP(-CI*CDT*TMP(J))
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      !$OMP PARALLEL DO PRIVATE(I,TMP,J)
      DO I = 1, NX
         RHO(I) = SUM(ABS(PHI(I,1:MS))**2)
         FZ(I)  = 2.0D0*ABS(PHI(I,1))**2 + ABS(PHI(I,2))**2&
                  - ABS(PHI(I,4))**2 - 2.0D0*ABS(PHI(I,5))**2
         TMP(1) = V(I) + TAU(0)*RHO(I) + 2.0D0*TAU(1)*FZ(I)&
                 + 0.4D0*TAU(2)*ABS(PHI(I,5))**2
         TMP(2) = V(I) + TAU(0)*RHO(I) + TAU(1)*FZ(I)&
                 + 0.4D0*TAU(2)*ABS(PHI(I,4))**2
         TMP(3) = V(I) + TAU(0)*RHO(I)+ 0.2D0*TAU(2)*ABS(PHI(I,3))**2
         TMP(4) = V(I) + TAU(0)*RHO(I) - TAU(1)*FZ(I) + 0.4D0*TAU(2)*ABS(PHI(I,2))**2 
         TMP(5) = V(I)+TAU(0)*RHO(I) - 2.0D0*TAU(1)*FZ(I)&
                  + 0.4D0*TAU(2)*ABS(PHI(I,1))**2
         DO J = 1, MS
            PHI(I,J) = PHI(I,J)*EXP(-CI*CDT*TMP(J))
         END DO
      END DO
      !$OMP END PARALLEL DO
  END SELECT  
END SUBROUTINE SP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                 Calculates FX, FY and FZ                                         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
SUBROUTINE FXYZ(PHI, FX, FY, FZ)
  USE BASIC_DATA, ONLY :  NX 
  USE CGPE_DATA, ONLY : SPIN, SX, SY, SZ
  USE DOUBLE_PRECISION
  !USE OMP_LIB

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
  REAL(KIND=DBL), DIMENSION(1:NX), INTENT(OUT) :: FX, FY, FZ
  INTEGER :: I, MS
  
  MS = 2*SPIN+1

  !$OMP PARALLEL DO PRIVATE(I)
  DO I = 1, NX
     FX(I) = REAL(DOT_PRODUCT(PHI(I,1:MS), MATMUL(SX(1:MS,1:MS),PHI(I,1:MS))))
     FY(I) = REAL(DOT_PRODUCT(PHI(I,1:MS), MATMUL(SY(1:MS,1:MS),PHI(I,1:MS))))
     FZ(I) = DOT_PRODUCT(SZ(1:MS),ABS(PHI(I,1:MS))**2)
  END DO
  !$OMP END PARALLEL DO 
END SUBROUTINE FXYZ

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            Calculates C12, C13, C23, C34, C35 and C45, i.e. the elements of Hamiltonian          !
!            corresponding to off-diagonal interaction terms in spin-2 systems                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
  USE BASIC_DATA, ONLY : NX
  USE CGPE_DATA, ONLY : TAU
  USE DOUBLE_PRECISION
  !USE OMP_LIB

  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:5), INTENT(IN) :: PHI
  COMPLEX(KIND=DBL), DIMENSION(1:NX), INTENT(IN) :: FMINUS
  COMPLEX(KIND=DBL), DIMENSION(1:NX), INTENT(OUT) :: C12, C13, C23, C34, C35, C45

  !$OMP PARALLEL WORKSHARE
  C12(:) = TAU(1)*FMINUS(:) - 0.4D0*TAU(2)*PHI(:,4)*CONJG(PHI(:,5))

  C13(:) = 0.2D0*TAU(2)*PHI(:,3)*CONJG(PHI(:,5))

  C23(:) = SQRT(1.5D0)*TAU(1)*FMINUS(:)-0.2D0*TAU(2)*PHI(:,3)*CONJG(PHI(:,4))
             
  C34(:) = SQRT(1.5D0)*TAU(1)*FMINUS(:)-0.2D0*TAU(2)*PHI(:,2)*CONJG(PHI(:,3)) 
            
  C35(:) = 0.2D0*TAU(2)*PHI(:,1)*CONJG(PHI(:,3))

  C45(:) = TAU(1)*FMINUS(:) - 0.4D0*TAU(2)*PHI(:,1)*CONJG(PHI(:,2))
  !$OMP END PARALLEL WORKSHARE
END SUBROUTINE MAT_C

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          Calculates the norm of individual components                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NORMC(PHI, NORM)
  USE BASIC_DATA, ONLY : NX, DX
  USE CGPE_DATA, ONLY : SPIN
  USE DOUBLE_PRECISION
  !USE OMP_LIB
  
  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: NORM

  INTERFACE
    PURE FUNCTION SIMPSON(F, DX)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: F
      REAL(KIND=DBL), INTENT(IN) :: DX
      REAL(KIND=DBL) :: SIMPSON
    END FUNCTION SIMPSON
  END INTERFACE
        
  REAL(KIND=DBL), DIMENSION(1:NX,1:2*SPIN+1) :: TMP1D
  INTEGER :: J, K, MS
  
  MS = 2*SPIN+1

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(J,K)
  DO K = 1, MS
     DO J = 1, NX
        TMP1D(J,K) = ABS(PHI(J,K))**2
     END DO
  END DO
  !$OMP END PARALLEL DO

  DO J = 1, MS
    NORM(J) = SIMPSON(TMP1D(:,J),DX)
  END DO
END SUBROUTINE NORMC

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Calculates the total norm and normalizes the total density to 1                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NORMT(PHI, TNORM)
  USE BASIC_DATA, ONLY : NX, DX
  USE CGPE_DATA, ONLY : SPIN
  USE DOUBLE_PRECISION
  !USE OMP_LIB

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(INOUT) :: PHI
  REAL(KIND=DBL), INTENT(OUT) :: TNORM

  INTERFACE
    PURE FUNCTION SIMPSON(F, DX) 
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: F
      REAL(KIND=DBL), INTENT(IN) :: DX
      REAL(KIND=DBL) :: SIMPSON
    END FUNCTION SIMPSON
  END INTERFACE

  INTEGER :: I, MS
  REAL(KIND=DBL), DIMENSION(1:NX) :: TMP1D

  MS = 2*SPIN+1

  !$OMP PARALLEL DO PRIVATE(I)
  DO I = 1, NX
     TMP1D(I) = SUM(ABS(PHI(I,1:MS))**2) 
  END DO
  !$OMP END PARALLEL DO
  TNORM = SQRT(SIMPSON(TMP1D, DX))
  PHI = PHI/TNORM
END SUBROUTINE NORMT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Calculates the root mean square sizes of the (2*SPIN+1) components                  !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE RAD(PHI, RMS)
  USE BASIC_DATA, ONLY : NX, DX
  USE CGPE_DATA, ONLY : X2, SPIN
  USE DOUBLE_PRECISION
  !USE OMP_LIB

  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: RMS

  INTERFACE
    PURE FUNCTION SIMPSON(F, DX)
    USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: F
      REAL(KIND=DBL), INTENT(IN) :: DX
      REAL(KIND=DBL) :: SIMPSON
    END FUNCTION SIMPSON
  END INTERFACE

  INTEGER :: I, J, MS
  REAL(KIND=DBL), DIMENSION(1:NX,1:2*SPIN+1) :: TMP1D

  MS = 2*SPIN+1

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J) 
  DO J = 1, MS
     DO I = 1, NX
        TMP1D(I,J) = X2(I)*ABS(PHI(I,J))**2
     END DO
  END DO
  !$OMP END PARALLEL DO

  DO J = 1, MS
     RMS(J) = SQRT(SIMPSON(TMP1D(:,J), DX))
  END DO
END SUBROUTINE RAD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       Calculates the (2*SPIN+1) component chemical potentials, energy, and magnetization         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE ENERGY(PHI, MU, EN, MZ)
  USE BASIC_DATA, ONLY : NX, DX, CI, EPS
  USE CGPE_DATA, ONLY : V, TAU, SPIN
  USE SOC_DATA , ONLY : GAMMAX
  USE DOUBLE_PRECISION
  !USE OMP_LIB

  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: MU
  REAL(KIND=DBL), INTENT(OUT) :: EN, MZ

  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: P
      REAL(KIND=DBL), INTENT(IN) :: DX
      REAL(KIND=DBL), DIMENSION(1:SIZE(P)) :: DP
    END FUNCTION DIFF

    PURE FUNCTION SIMPSON(F, DX)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: F
      REAL(KIND=DBL), INTENT(IN) :: DX
      REAL(KIND=DBL) :: SIMPSON
    END FUNCTION SIMPSON

    SUBROUTINE FXYZ(PHI, FX, FY, FZ)
      USE BASIC_DATA, ONLY : NX, CI
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(1:NX), INTENT(OUT) :: FX, FY, FZ
    END SUBROUTINE FXYZ

    SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      USE BASIC_DATA, ONLY : NX
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:5), INTENT(IN) :: PHI
      COMPLEX(KIND=DBL), DIMENSION(1:NX), INTENT(IN) ::  FMINUS
      COMPLEX(KIND=DBL), DIMENSION(1:NX), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
    END SUBROUTINE MAT_C

    SUBROUTINE NORMC(PHI, NORM)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: NORM
    END  SUBROUTINE NORMC
  END INTERFACE
!---------------------------------------------------------------------------------------------------

  INTEGER :: K, I
  REAL(KIND=DBL), DIMENSION(1:NX) ::  FX, FY, FZ, RHO
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:2*SPIN+1) ::  DPSX
  REAL(KIND=DBL), DIMENSION(1:NX,1:2*SPIN+1) :: DPX
  COMPLEX(KIND=DBL), DIMENSION(1:NX) :: THETA, FMINUS, C12, C13, C23, C34, C35, C45
  REAL(KIND=DBL), DIMENSION(1:2*SPIN+1) :: NORM
  REAL(KIND=DBL), DIMENSION(1:NX, 1:2*SPIN+1+2) :: TMP1D

  CALL FXYZ(PHI, FX, FY, FZ)

  CALL NORMC(PHI, NORM)

  DO K = 1, 2*SPIN+1
     DPSX(1:NX,K) = (DIFF(REAL(PHI(1:NX,K)), DX) + CI*DIFF(AIMAG(PHI(1:NX,K)), DX))
  END DO

  DPX(:,:) = ABS(DPSX(:,:))**2
 
  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO PRIVATE(I)
      DO I = 1, NX
         RHO(I) = SUM(ABS(PHI(I,1:2*SPIN+1))**2) 
         TMP1D(I,1) = REAL(V(I)*RHO(I) + TAU(0)*RHO(I)**2/2.0D0 +&
                      TAU(1)*((FX(I))**2 +(FY(I))**2 +(FZ(I))**2)/2.0D0+&
                      DPX(I,1) + DPX(I,2) + DPX(I,3) +&
                      SQRT(2.0D0)*GAMMAX*(-CI*CONJG(PHI(I,1))*DPSX(I,2)-&
                      CI*CONJG(PHI(I,3))*DPSX(I,2)- &
                      CI*CONJG(PHI(I,2))*(DPSX(I,1)+DPSX(I,3))))

         TMP1D(I,2) = ABS(PHI(I,1))**2-ABS(PHI(I,3))**2

         TMP1D(I,3) = REAL(V(I)*ABS(PHI(I,1))*ABS(PHI(I,1)) + DPX(I,1) +&
                      TAU(0)*(ABS(PHI(I,1))**2+ABS(PHI(I,2))**2+ABS(PHI(I,3))**2)*ABS(PHI(I,1))**2+&
                      TAU(1)*(ABS(PHI(I,1))**2+ABS(PHI(I,2))**2-ABS(PHI(I,3))**2)*ABS(PHI(I,1))**2+&
                      TAU(1)*CONJG(PHI(I,3))*PHI(I,2)**2*CONJG(PHI(I,1))+&
                      SQRT(2.0D0)*GAMMAX*(-CI*CONJG(PHI(I,1))*DPSX(I,2)))

         TMP1D(I,4) = REAL(V(I)*ABS(PHI(I,2))*ABS(PHI(I,2)) + DPX(I,2) +&
                      TAU(0)*(ABS(PHI(I,1))**2+ABS(PHI(I,2))**2+ABS(PHI(I,3))**2)*ABS(PHI(I,2))**2+&
                      TAU(1)*(ABS(PHI(I,1))**2+ABS(PHI(I,3))**2)*ABS(PHI(I,2))**2+&
                      2.0D0*TAU(1)*PHI(I,1)*PHI(I,3)*CONJG(PHI(I,2))**2+&
                      SQRT(2.0D0)*GAMMAX*(-CI*CONJG(PHI(I,2))*(DPSX(I,1)+DPSX(I,3))))

         TMP1D(I,5) = REAL(V(I)*ABS(PHI(I,3))*ABS(PHI(I,3)) + DPX(I,3) +&
                      TAU(0)*(ABS(PHI(I,1))**2+ABS(PHI(I,2))**2+ABS(PHI(I,3))**2)*ABS(PHI(I,3))**2+&
                      TAU(1)*(ABS(PHI(I,2))**2+ABS(PHI(I,3))**2-ABS(PHI(I,1))**2)*ABS(PHI(I,3))**2+&
                      TAU(1)*CONJG(PHI(I,1))*PHI(I,2)**2*CONJG(PHI(I,3))+&
                      SQRT(2.0D0)*GAMMAX*(-CI*CONJG(PHI(I,3))*DPSX(I,2)))

      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      FMINUS = FX - CI * FY
      CALL MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      !$OMP PARALLEL DO PRIVATE(I)
      DO I = 1, NX
         RHO(I) =  SUM(ABS(PHI(I,1:2*SPIN+1))**2)
         THETA(I) = (2.0D0*PHI(I,1)*PHI(I,5)-2.0D0*PHI(I,2)*PHI(I,4)+PHI(I,3)**2)/SQRT(5.0D0)
         TMP1D(I,1) = REAL(V(I)*RHO(I) + TAU(0)*RHO(I)**2/2.0D0 +&
                      TAU(1)*((FX(I))**2 +(FY(I))**2 +(FZ(I))**2)/2.0D0+&
                      TAU(2)*(THETA(I)*CONJG(THETA(I)))/2.0D0+&
                      DPX(I,1) + DPX(I,2) + DPX(I,3) + DPX(I,4) + DPX(I,5)+&
                      GAMMAX*(-CI*2.0D0*(CONJG(PHI(I,1))*DPSX(I,2)) - CI*2.0D0*(CONJG(PHI(I,2))*DPSX(I,1)) &
                      -CI*SQRT(6.0D0)*(CONJG(PHI(I,2))*DPSX(I,3) + CONJG(PHI(I,3))*DPSX(I,2) +&
                      CONJG(PHI(I,3))*DPSX(I,4)+ CONJG(PHI(I,4))*DPSX(I,3))-&
                      CI*2.0D0*(CONJG(PHI(I,4))*DPSX(I,5) + CONJG(PHI(I,5))*DPSX(I,4))))

         TMP1D(I,2) = 2.0D0*ABS(PHI(I,1))**2 + ABS(PHI(I,2))**2 - ABS(PHI(I,4))**2-&
                      2.0D0*ABS(PHI(I,5))**2

         TMP1D(I,3) = REAL((V(I)+TAU(0)*RHO(I) + 2.0D0*TAU(1)*FZ(I)+&
                      0.4D0*TAU(2)*ABS(PHI(I,5))**2)*ABS(PHI(I,1))**2+&
                      (C12(I)*PHI(I,2)+C13(I)*PHI(I,3))*CONJG(PHI(I,1))&
                      -CI*2.0D0*GAMMAX*(CONJG(PHI(I,1))*DPSX(I,2))+&
                      ABS(DPSX(I,1))**2)

         TMP1D(I,4) = REAL((V(I)+TAU(0)*RHO(I) + TAU(1)*FZ(I)+&
                      0.4D0*TAU(2)*ABS(PHI(I,4))**2)*ABS(PHI(I,2))**2+&
                      (CONJG(C12(I))*PHI(I,1) + C23(I)*PHI(I,3))*CONJG(PHI(I,2))+&
                      GAMMAX*(-CI*2.0D0*(CONJG(PHI(I,2))*DPSX(I,1))-CI*SQRT(6.0D0)*CONJG(PHI(I,2))*DPSX(I,3))&
                      +ABS(DPSX(I,2))**2)

         TMP1D(I,5) = REAL((V(I)+TAU(0)*RHO(I) +&
                      0.2D0*TAU(2)*ABS(PHI(I,3))**2)*ABS(PHI(I,3))**2+&
                      (CONJG(C13(I))*PHI(I,1)+CONJG(C23(I))*PHI(I,2) + C34(I)*PHI(I,4)+C35(I)*PHI(I,5))*CONJG(PHI(I,3))&
                      -CI*GAMMAX*SQRT(6.0D0)*(CONJG(PHI(I,3))*DPSX(I,2) + CONJG(PHI(I,3))*DPSX(I,4)) &
                      +ABS(DPSX(I,3))**2)

         TMP1D(I,6) = REAL((V(I)+TAU(0)*RHO(I) - TAU(1)*FZ(I)+&
                      0.4D0*TAU(2)*ABS(PHI(I,2))**2)*ABS(PHI(I,4))**2+&
                      (CONJG(C34(I))*PHI(I,3)+C45(I)*PHI(I,5))*CONJG(PHI(I,4))+&
                      GAMMAX*(-CI*SQRT(6.0D0)*CONJG(PHI(I,4))*DPSX(I,3) -CI*2.0D0*CONJG(PHI(I,4))*DPSX(I,5))+&
                      ABS(DPSX(I,4))**2)

         TMP1D(I,7) = REAL((V(I)+TAU(0)*RHO(I) - 2.0D0*TAU(1)*FZ(I)+&
                      0.4D0*TAU(2)*ABS(PHI(I,1))**2)*ABS(PHI(I,5))**2+&
                      (CONJG(C35(I))*PHI(I,3)+CONJG(C45(I))*PHI(I,4))*CONJG(PHI(I,5))&
                      -CI*2.0D0*GAMMAX*CONJG(PHI(I,5))*DPSX(I,4)+&
                      ABS(DPSX(I,5))**2)
      END DO
      !$OMP END PARALLEL DO
  END SELECT
  EN = SIMPSON(TMP1D(:,1),DX)/SUM(NORM)  
  MZ = SIMPSON(TMP1D(:,2),DX)
  DO K = 1, 2*SPIN+1
     MU(K) = SIMPSON(TMP1D(:,K+2),DX)/(NORM(K)+ EPS)
  END DO
END SUBROUTINE ENERGY 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                  Solves the non-linear equations by using Newton Raphson method                  !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NEWTON(X, N1, N2, N3, N4, N5, MAG, ERR_MSG)
  USE BASIC_DATA, ONLY : TOL_NR
  USE CGPE_DATA, ONLY : NTRIAL
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  REAL(KIND=DBL), INTENT(IN) :: N1, N2, N3, N4, N5, MAG
  INTEGER, INTENT(INOUT) :: ERR_MSG
  REAL(KIND=DBL), DIMENSION(1:2), INTENT(INOUT) :: X

  !USES USRFUN
  !Given an initial guess x for a root in n dimensions, take ntrial Newton-Raphson steps to
  !improve the root. 

  INTEGER :: I, K
  REAL(KIND=DBL) :: FJAC(2,2), FVEC(2), P(2)

  INTERFACE
    SUBROUTINE USRFUN(X,FVEC,FJAC,N1,N2,N3,N4,N5,MAG)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:2), INTENT(IN) :: X
      REAL(KIND=DBL), INTENT(OUT) :: FJAC(1:2,1:2)
      REAL(KIND=DBL), INTENT(OUT) :: FVEC(1:2)
      REAL(KIND=DBL), INTENT(IN) :: N1, N2, N3, N4, N5, MAG
    END SUBROUTINE USRFUN
  END INTERFACE

  DO K = 1, NTRIAL
     CALL USRFUN(X,FVEC,FJAC,N1,N2,N3,N4,N5,MAG) !User subroutine supplies function values at x in fvec
                                                  !and Jacobian matrix in fjac. 

     P(2) = (FVEC(2)*FJAC(1,1)-FVEC(1)*FJAC(2,1))/(FJAC(2,1)*FJAC(1,2)-FJAC(1,1)*FJAC(2,2))
     P(1) = (-FVEC(1)-FJAC(1,2)*P(2))/FJAC(1,1)

     DO I = 1, 2                    !Update solution.
        X(I) = X(I) + P(I)
     END DO
     IF(MAXVAL(ABS(P)).LE.TOL_NR) EXIT !Check root convergence.
  END DO

  IF(MAXVAL(ABS(P)).GT.TOL_NR)THEN
     WRITE(100,*) 'Newton-Raphason subroutine could not calculate the roots within' , TOL_NR, ' tolerance'
     ERR_MSG = 1
  END IF
END
SUBROUTINE USRFUN(X,FVEC,FJAC,N1,N2,N3,N4,N5,MAG)
  USE BASIC_DATA, ONLY : DT
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  REAL(KIND=DBL), DIMENSION(1:2), INTENT(IN) :: X
  REAL(KIND=DBL), INTENT(OUT) :: FJAC(1:2,1:2)
  REAL(KIND=DBL), INTENT(OUT) :: FVEC(1:2)
  REAL(KIND=DBL), INTENT(IN) :: N1, N2, N3, N4, N5, MAG

  FVEC(1) = EXP(2.0D0*(X(1)+2.0D0*X(2))*DT)*N1 + EXP(2.0D0*(X(1)+X(2))*DT)*N2 +&
            EXP(2.0D0*X(1)*DT)*N3 + EXP(2.0D0*(X(1)-X(2))*DT)*N4 +&
            EXP(2.0D0*(X(1)-2.0D0*X(2))*DT)*N5-1.0D0
  FVEC(2) = 2.0D0*EXP(2.0D0*(X(1)+2.0D0*X(2))*DT)*N1 +  EXP(2.0D0*(X(1)+X(2))*DT)*N2-&
            EXP(2.0D0*(X(1)-X(2))*DT)*N4 - 2.0D0*EXP(2.0D0*(X(1)-2.0D0*X(2))*DT)*N5-MAG

  FJAC(1,1) = (EXP(2.0D0*(X(1)+2.0D0*X(2))*DT)*N1 + EXP(2.0D0*(X(1)+X(2))*DT)*N2 +&
            EXP(2.0D0*X(1)*DT)*N3 + EXP(2.0D0*(X(1)-X(2))*DT)*N4 +&
            EXP(2.0D0*(X(1)-2.0D0*X(2))*DT)*N5)*2.0D0*DT
  FJAC(1,2) = EXP(2.0D0*(X(1)+2.0D0*X(2))*DT)*N1*4.0D0*DT + EXP(2.0D0*(X(1)+X(2))*DT)*N2*2.0D0*DT &
            + EXP(2.0D0*(X(1)-X(2))*DT)*N4*(-2.0D0*DT) +&
            EXP(2.0D0*(X(1)-2.0D0*X(2))*DT)*N5*(-4.0D0*DT)
  FJAC(2,1) = (2.0D0*EXP(2.0D0*(X(1)+2.0D0*X(2))*DT)*N1 +  EXP(2.0D0*(X(1)+X(2))*DT)*N2-&
            EXP(2.0D0*(X(1)-X(2))*DT)*N4 - 2.0D0*EXP(2.0D0*(X(1)-2.0D0*X(2))*DT)*N5)*2.0D0*DT
  FJAC(2,2) = 2.0D0*EXP(2.0D0*(X(1)+2.0D0*X(2))*DT)*N1*4.0D0*DT +  EXP(2.0D0*(X(1)+X(2))*DT)*N2*2.0D0*DT-&
            EXP(2.0D0*(X(1)-X(2))*DT)*N4*(-2.0D0*DT) - 2.0D0*EXP(2.0D0*(X(1)-2.0D0*X(2))*DT)*N5*(-4.0D0*DT)
END

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                     Calculates the discrete forward Fourier transform                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE FFT()
  USE CGPE_DATA, ONLY : PHI, PHIF, SPIN
  USE FFTW3, ONLY : FFTW_EXECUTE_DFT
  USE FFTW_DATA

  IMPLICIT NONE
  INTEGER :: L, MS

  MS = 2*SPIN+1
  DO L = 1, MS
     FFTFX(:) = PHI(:,L)
     CALL FFTW_EXECUTE_DFT(PLANFX,FFTFX,FFTBX)
     PHIF(:,L) = FFTBX(:)
  END DO
END SUBROUTINE FFT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         Calculates the discrete backward Fourier transform of component wavefunctions            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE BFT()
  USE BASIC_DATA, ONLY : NX
  USE CGPE_DATA, ONLY : PHI, PHIF, SPIN
  USE FFTW3, ONLY : FFTW_EXECUTE_DFT
  USE FFTW_DATA

  IMPLICIT NONE
  INTEGER :: L, MS

  MS = 2*SPIN+1
  DO L = 1, MS
     FFTBX(:) = PHIF(:,L)
     CALL FFTW_EXECUTE_DFT(PLANBX,FFTBX,FFTFX)
     PHI(:,L) = FFTFX(:)/DBLE(NX)
  END DO
END SUBROUTINE BFT
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               Solves the Hamiltonian corresponding to off-diagonal interaction terms             !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SE()
  USE BASIC_DATA, ONLY : CI, NX, CDT, EPS
  USE CGPE_DATA , ONLY : PHI, SPIN,  TAU 
  USE DOUBLE_PRECISION
 
  IMPLICIT NONE 

  INTERFACE
    SUBROUTINE FXYZ(PHI, FX, FY, FZ)
      USE BASIC_DATA, ONLY : NX, CI
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(1:NX), INTENT(OUT) :: FX, FY, FZ
    END SUBROUTINE FXYZ

    SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      USE BASIC_DATA, ONLY : NX
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:5), INTENT(IN) :: PHI
      COMPLEX(KIND=DBL), DIMENSION(1:NX), INTENT(IN) :: FMINUS
      COMPLEX(KIND=DBL), DIMENSION(1:NX), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
    END SUBROUTINE MAT_C
  END INTERFACE
 
  INTEGER  :: I, K
  COMPLEX(KIND=DBL) :: A1, B1, OMEGA, DELTA, KAPPA
  COMPLEX(KIND=DBL), DIMENSION(1:3) :: TMP

  REAL(KIND=DBL), DIMENSION(1:NX) :: FX, FY, FZ, FX1, FY1, FZ1
  COMPLEX(KIND=DBL), DIMENSION(1:NX) :: FMINUS, C12, C13, C23, C34, C35, C45,&
                                 C121, C131, C231, C341, C351, C451, FMINUS1
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:5) :: PHI1

  EXTERNAL :: ZHEEV
  INTEGER, PARAMETER :: LDA = 5
  INTEGER, PARAMETER :: LWMAX = 1000
  INTEGER :: INFO , LWORK
  REAL(KIND=DBL) :: RWORK(13), W(5) 
  COMPLEX(KIND=DBL) :: A(1:5, 1:5), WORK(LWMAX)
    
  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO PRIVATE(I,TMP,A1,B1,OMEGA,DELTA,KAPPA) 
      DO I = 1, NX
         A1 =  TAU(1) * PHI(I,2) * CONJG(PHI(I,3))
         B1 =  TAU(1) * PHI(I,2) * CONJG(PHI(I,1))
         OMEGA = CDT*SQRT((ABS(A1))**2+(ABS(B1))**2)

         IF(ABS(OMEGA).GE.EPS)THEN
            DELTA = (COS(OMEGA) -1)/(OMEGA)**2
            KAPPA = -CI * SIN(OMEGA)/OMEGA
         ELSE
            DELTA = -0.5D0
            KAPPA = -CI
         END IF

         TMP(1) = (DELTA*CDT**2)*(A1*CONJG(A1)*PHI(I,1)+&
                  A1*CONJG(B1)*PHI(I,3)) + KAPPA*CDT*A1&
                  *PHI(I,2)
         TMP(2) = (DELTA*CDT**2)*(A1*CONJG(A1)+B1*CONJG(B1))&
                  *PHI(I,2) + KAPPA*CDT*&
                  (CONJG(A1)*PHI(I,1)+CONJG(B1)*PHI(I,3))
         TMP(3) = (DELTA*CDT**2)*(CONJG(A1)*B1*&
                  PHI(I,1)+B1*CONJG(B1)*PHI(I,3)) + KAPPA*CDT*&
                  B1*PHI(I,2)

         PHI(I,1) = PHI(I,1) + TMP(1)
         PHI(I,2) = PHI(I,2) + TMP(2)
         PHI(I,3) = PHI(I,3) + TMP(3)
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      CALL FXYZ(PHI, FX, FY, FZ)

      FMINUS(:) = FX(:) - CI * FY(:)

      CALL MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)

      PHI1(:,1) = PHI(:,1) - (CI*CDT/2.0D0)*(C12(:)*PHI(:,2)+C13(:)*PHI(:,3))
      PHI1(:,2) = PHI(:,2) - (CI*CDT/2.0D0)*(CONJG(C12(:))*PHI(:,1)+C23(:)*PHI(:,3))
      PHI1(:,3) = PHI(:,3) - (CI*CDT/2.0D0)*(CONJG(C13(:))*PHI(:,1)+CONJG(C23(:))*PHI(:,2)+C34(:)*PHI(:,4)+C35(:)*PHI(:,5))
      PHI1(:,4) = PHI(:,4) - (CI*CDT/2.0D0)*(CONJG(C34(:))*PHI(:,3)+C45(:)*PHI(:,5))
      PHI1(:,5) = PHI(:,5) - (CI*CDT/2.0D0)*(CONJG(C35(:))*PHI(:,3)+CONJG(C45(:))*PHI(:,4))

      CALL FXYZ(PHI1, FX1, FY1, FZ1)

      FMINUS1 = FX1 - CI * FY1

      CALL MAT_C(PHI1, FMINUS1, C121, C131, C231, C341, C351, C451)

      C12 = (C12+C121)/2.0D0
      C13 = (C13+C131)/2.0D0
      C23 = (C23+C231)/2.0D0
      C34 = (C34+C341)/2.0D0
      C35 = (C35+C351)/2.0D0
      C45 = (C45+C451)/2.0D0

      A = CMPLX(0.0D0,0.0D0,KIND=DBL)

      !$OMP PARALLEL DO PRIVATE(I,W,WORK,LWORK,RWORK,INFO) FIRSTPRIVATE(A)
      DO I = 1, NX
         A(1,2) = C12(I)
         A(1,3) = C13(I)
    
         A(2,1) = CONJG(C12(I))
         A(2,3) = C23(I)

         A(3,1) = CONJG(C13(I))
         A(3,2) = CONJG(C23(I))
         A(3,4) = C34(I)
         A(3,5) = C35(I)
 
         A(4,3) = CONJG(C34(I))
         A(4,5) = C45(I)

         A(5,3) = CONJG(C35(I))
         A(5,4) = CONJG(C45(I))

         LWORK = -1
         CALL ZHEEV( 'Vectors', 'Lower', 5, A, LDA, W, WORK, LWORK, RWORK, INFO )
         LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
         CALL ZHEEV( 'Vectors', 'Lower', 5, A, LDA, W, WORK, LWORK, RWORK, INFO )

         PHI(I,1:5) = MATMUL(CONJG(TRANSPOSE(A(1:5,1:5))), PHI(I,1:5))
         DO K = 1, 5
            PHI(I,K) = EXP(-CI*CDT*W(K)) * PHI(I,K)
         END DO
         PHI(I,1:5) = MATMUL(A(1:5,1:5), PHI(I,1:5))   
   
         A = CMPLX(0.0D0,0.0D0,KIND=DBL) 
         W = 0.0D0 
      END DO
      !$OMP END PARALLEL DO 
  END SELECT
END SUBROUTINE SE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       Solves the Hamiltonian corresponding to spin-orbit coupling terms in Fourier space         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SOC()
  USE BASIC_DATA, ONLY : CI, NX, CDT
  USE CGPE_DATA, ONLY : PHIF, KX, SPIN
  USE SOC_DATA, ONLY : GAMMAX, A
  !USE OMP_LIB

  IMPLICIT NONE
  INTEGER :: I, J, MS

  MS = 2*SPIN+1
  !$OMP PARALLEL DO PRIVATE(I,J)
  DO I = 1, NX
     PHIF(I,1:MS) = MATMUL(TRANSPOSE(A(1:MS,1:MS)), PHIF(I,1:MS))
     DO J = -SPIN, SPIN, 1
        PHIF(I,J+SPIN+1) = EXP(-CI*CDT*(2.0D0*J*GAMMAX*KX(I))) * PHIF(I,J+SPIN+1)
     END DO
     PHIF(I,1:MS) = MATMUL(A(1:MS,1:MS), PHIF(I,1:MS))
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE SOC

