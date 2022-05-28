! File name: cgpe2d.f90 
! Title: FORTRESS: FORTRAN programs to solve coupled
! Gross-Pitaevskii equations for spin-orbit coupled spin-f
! Bose-Einstein condensate with spin f = 1 or 2
!        
! Authors: Paramjeet Banger, Pardeep Kaur, Arko Roy, & Sandeep Gautam
!
! Fortran code to solve 2D coupled Gross-Pitaevskii equations (CGPEs)
! for spin-2(1) BEC using time-splitting spectral method.       
!
!---------------------------------------------------------------------------------------------------
! PI = 4.0*arctan(1), CI = \sqrt{-1}                                                          
! NITER = Number of total iterations 
! NSTP  = Number of iterations after which transeint results are written
! OPENMP_THREADS = Number of threads chosen for OpenMP parallelization
! FFTW_THREADS = Number of FFTW threads
! NX = Number of space grid points in x-direction
! NY = Number of space grid points in y-direction
! DX = Spatial step size along x
! DY = Spatial step size along y
! DT = Temporal step size
! LX = NX*DX, LY = NY*DY, (2D spatial domain)
! AMU = atomic mass unit in Kg
! HBAR = reduced Planck's constat in SI units
! CDT = Complex varible which will be defined as '-idt' and 'dt' in
!       imaginary and realtime propagations, respectively 
!---------------------------------------------------------------------------------------------------

MODULE DOUBLE_PRECISION
  INTEGER, PARAMETER :: DBL = KIND(0.0D0)
END MODULE DOUBLE_PRECISION

MODULE BASIC_DATA
  USE DOUBLE_PRECISION
  REAL(KIND=DBL), PARAMETER :: PI = 4.0D0*DATAN(1.0D0)
  COMPLEX(KIND=DBL), PARAMETER :: CI = CMPLX(0.0D0,1.0D0,KIND=DBL)
  REAL(KIND=DBL), PARAMETER :: AMU = 1.66D-27, HBAR = 1.054560653D-34, AU = 0.529177208D-10
  REAL(KIND=DBL), PARAMETER :: EPS = 1.0D-40
  REAL(KIND=DBL), PARAMETER :: TOL_NR = 1.0D-6
  INTEGER :: NITER
  INTEGER :: OPENMP_THREADS, FFTW_THREADS
  INTEGER :: NX, NX2, NXX, NY, NY2, NYY
  REAL(KIND=DBL) :: DX, DY, DT
  REAL(KIND=DBL) :: LX, LY
  INTEGER :: STP, NSTP
  COMPLEX(KIND=DBL) :: CDT
END MODULE BASIC_DATA

MODULE CGPE_DATA
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : NX, NY
!---------------------------------------------------------------------------------------------------
  !Scattering length values for Na-23-Antiferromagntic
  !REAL(KIND=DBL), PARAMETER :: M = (23.0D0*AMU)
  !REAL(KIND=DBL), PARAMETER :: A0 = 34.9D0*AU, A2 = 45.8D0*AU, A4 = 64.5D0*AU
  !Scattering length values for Rb-83-Ferromagnetic
  !REAL(KIND=DBL), PARAMETER :: M = (83.0D0*AMU) 
  !REAL(KIND=DBL), PARAMETER :: A0 = 83.0D0*AU, A2 = 82.0D0*AU, A4 = 81.0D0*AU
  !Scattering length values for Rb-87-Cyclic
  !REAL(KIND=DBL), PARAMETER :: M = (87.0D0*AMU) 
  !REAL(KIND=DBL), PARAMETER :: A0 = 87.93D0*AU, A2 = 91.28D0*AU, A4 = (99.18D0)*AU 
!---------------------------------------------------------------------------------------------------
  REAL(KIND=DBL) :: M
  REAL(KIND=DBL) :: A0, A2, A4
  REAL(KIND=DBL) :: NUX, NUY, NUZ,&
                    ALPHAX, ALPHAY, ALPHAZ
  INTEGER :: NATOMS, SPIN 
  INTEGER, PARAMETER :: NTRIAL = 1000
  REAL(KIND=DBL) :: TOL
  INTEGER :: SWITCH_IM, OPTION_FPC
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: X, X2
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: Y, Y2
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: KX
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: KY
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:) :: V, R2
  COMPLEX(KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:) :: PHI, PHIF, PHI_OLD
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
  REAL(KIND=DBL) :: GAMMAX,  GAMMAY
  INTEGER :: SWITCH_SOC
END MODULE SOC_DATA

MODULE FFTW_DATA
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : NX, NY
  USE FFTW3, ONLY : C_DOUBLE_COMPLEX, C_PTR, C_INT
  COMPLEX(KIND=C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:) :: FFTFXY, FFTBXY
  TYPE(C_PTR) :: PLANFX, PLANBX, PLANFXY, PLANBXY
  INTEGER(KIND=C_INT) :: THREADS_INIT
END MODULE FFTW_DATA
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                        Main Program - CGPE2D                                     ! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PROGRAM CGPE2D
  USE BASIC_DATA
  USE CGPE_DATA
  USE SOC_DATA
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
      COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: NORM
    END SUBROUTINE NORMC
 
    SUBROUTINE NORMT(PHI, TNORM)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(INOUT) :: PHI
      REAL(KIND=DBL), INTENT(OUT) :: TNORM
    END SUBROUTINE NORMT

    SUBROUTINE RAD(PHI, RMS)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: RMS
    END SUBROUTINE RAD

    SUBROUTINE ENERGY(PHI, MU, EN, MZ)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
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
  
  INTEGER :: I, J, K, MS, IOSTAT, II, ERR_MSG
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

  INQUIRE (FILE='./INPUT/input2D.dat', EXIST = FEXIST)
  IF(FEXIST)THEN
     OPEN (10 ,FILE='./INPUT/input2D.dat', STATUS='OLD',FORM='FORMATTED',&
           ACTION='READ')
  ELSE
     WRITE(100,*) "input2D.dat' does not exist."
     PRINT*, "input2D.dat' does not exist."
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
  READ(10,*) NX, NY
  READ(10,*) DX, DY, DT
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
  READ(10,*) GAMMAX, GAMMAY
  IF(ABS(GAMMAX).GT.EPS .OR. ABS(GAMMAY).GT.EPS) SWITCH_SOC = 1
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
  NY2 = NY/2+1
  NYY = NY-1

  ALPHAX = NUX/NUX
  ALPHAY = NUY/NUX
  ALPHAZ = NUZ/NUX

  LX = NX*DX
  LY = NY*DY
  
  ALLOCATE(NORM(1:2*SPIN+1), RMS(1:2*SPIN+1), MU(1:2*SPIN+1), SIGMA(1:2*SPIN+1))
  ALLOCATE(TMP(1:4*SPIN+2))
  CALL ALLOCATE_MEM()

  !$ CALL OMP_set_num_threads(OPENMP_THREADS)
  !$OMP PARALLEL
    !$OMP SINGLE
    !$ WRITE(*,*) 'MAXIMUM_THREAD_NUMBER_SET = ', OMP_get_num_threads()
    !$OMP END SINGLE
  !$OMP END PARALLEL

  IF (SPIN.EQ.1)THEN
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
      OPEN (9, FILE='./OUTPUT/convergence.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
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
      TAU(0) = 2.0D0*SQRT(ALPHAZ/(2.0D0*PI))*(4.0D0*PI*REAL(NATOMS)*(A0*AU+2.0D0*A2*AU)/(3.0D0*AOSC))
      TAU(1) = 2.0D0*SQRT(ALPHAZ/(2.0D0*PI))*(4.0D0*PI*REAL(NATOMS)*(A2*AU-A0*AU)/(3.0D0*AOSC))
    CASE(2)
      TAU(0) = 2.0D0*(SQRT(ALPHAZ/(2.0D0*PI))*4.0D0*PI*DBLE(NATOMS)*(4.0D0*A2*AU+3.0D0*A4*AU)/(7.0D0*AOSC))
      TAU(1) = 2.0D0*(SQRT(ALPHAZ/(2.0D0*PI))*4.0D0*PI*DBLE(NATOMS)*(A4*AU-A2*AU)/(7.0D0*AOSC))
      TAU(2) = 2.0D0*(SQRT(ALPHAZ/(2.0D0*PI))*4.0D0*PI*DBLE(NATOMS)*(7.0D0*A0*AU-10.0D0*A2*AU+3.0D0*A4*AU)/(7.0D0*AOSC))
  END SELECT

  WRITE(2,*)
  WRITE(2,899) OPENMP_THREADS, FFTW_THREADS
  WRITE(2,900) SWITCH_IM, OPTION_FPC, SWITCH_SOC, GAMMAX, GAMMAY, TOL
  WRITE(2,901) ALPHAX, ALPHAY, ALPHAZ
  WRITE(2,*)
  WRITE(2,902) NX, NY
  WRITE(2,903) NITER, NSTP, STP
  WRITE(2,905) DX,  DY
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

  899 FORMAT(' OPENMP_THREADS = ', I3,', FFTW_THREADS = ', I3)
  900 FORMAT(' SWITCH_IM = ', I2,', OPTION_FPC = ', I2, ', SWITCH_SOC = ',&
                I2, ', GAMMAX = ', F6.2, ', GAMMAY = ', F6.2, ', TOL = ', ES9.2)
  901 FORMAT(' Anisotropy ALPHAX =',F12.6,', ALPHAY = ',F12.6, ', ALPHAZ = ',F12.6)
  902 FORMAT(' No. of space steps NX = ',I8, ', NY = ', I8)
  903 FORMAT(' No. of time steps : NITER = ',I9,', NSTP = ',I9,', STP = ', I9)
  904 FORMAT(' A0 = ', F12.6, ', A2 = ', F12.6, ',  (all in units of Bohr radius) ')
  905 FORMAT(' Space Step DX = ',F10.6,', DY = ',F10.6)
  906 FORMAT(' Time Step  DT = ', F10.6,', MAG = ',F12.6)
  907 FORMAT(' AOSC = ', E12.4 )
  908 FORMAT(' TAU(0) = ', F12.6,', TAU(1) = ', F12.6,' (in dimensionless units)')
  911 FORMAT(' A0 = ', F12.6, ', A2 = ', F12.6, ', A4 = ', F12.6, ' (all in units of Bohr radius) ')
  912 FORMAT(' TAU(0) = ', F12.6,', TAU(1) = ', F12.6,', TAU(2) = ', F12.6, ' (in dimensionless units)')

  CALL INITIALIZE()
    
  IF((SWITCH_IM.NE.1))THEN
    IF(IOSTAT.NE.0) THEN
      PRINT*, 'Check if the "initial_sol_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" is available'
      PRINT*,'if not, user can read the "solution_file_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" file as an initial solution.'
      WRITE(100,*) 'Check if the "initial_sol_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" is available' 
      WRITE(100,*) 'if not, user can read the "solution_file_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat" file as an initial solution.'
      STOP
    END IF       
    DO I = 1, NX
       DO J = 1, NY
          READ(1,*) X(I), Y(J), (TMP(K), K = 1, 4*SPIN+2)
          SELECT CASE(SPIN)
            CASE(1)
              PHI(I,J,1) = SQRT(TMP(1))*EXP(CI*TMP(2))
              PHI(I,J,2) = SQRT(TMP(3))*EXP(CI*TMP(4))
              PHI(I,J,3) = SQRT(TMP(5))*EXP(CI*TMP(6))
            CASE(2)       
              PHI(I,J,1) = SQRT(TMP(1))*EXP(CI*TMP(2))
              PHI(I,J,2) = SQRT(TMP(3))*EXP(CI*TMP(4))
              PHI(I,J,3) = SQRT(TMP(5))*EXP(CI*TMP(6))
              PHI(I,J,4) = SQRT(TMP(7))*EXP(CI*TMP(8))
              PHI(I,J,5) = SQRT(TMP(9))*EXP(CI*TMP(10))
          END SELECT
       END DO
    END DO
  END IF
 
  IF (SWITCH_IM .EQ. 1) CALL NORMT(PHI, TNORM)

  MS = 2*SPIN+1

  CALL NORMC(PHI, NORM)
  CALL RAD(PHI, RMS)
  CALL ENERGY(PHI, MU, EN, MZ)
  
  SELECT CASE(SPIN)
    CASE(1)
      WRITE (2, 1001)
      WRITE (2, 1002)
      WRITE (2, 1001)
      WRITE (2, 1003) SUM(NORM), EN/2.0D0, (MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,NY2,K)), K = 1, MS)
    CASE(2)
      WRITE (2, 1011)
      WRITE (2, 1012)
      WRITE (2, 1011)
      WRITE (2, 1013) SUM(NORM), EN/2.0D0, (MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,NY2,K)), K = 1, MS)
  END SELECT

  1001 FORMAT (23X, 76('-'))
  1002 FORMAT (23X, 'Norm', 7X,'EN', 7X, 'MU1', 7X, 'MU2', 7X, 'MU3', 6X,&
                    'phi1', 6X, 'phi2', 6X,'phi3')
  1003 FORMAT ('Initial: ', 8X, 1F11.4, 4F10.5, 3F10.5)
  1011 FORMAT (22X, 117('-'))
  1012 FORMAT (23X, 'Norm', 7X,'EN', 7X, 'MU1', 7X, 'MU2', 7X, 'MU3', 7X, 'MU4', 7X, 'MU5', 7X,&
                    'phi1', 7X, 'phi2', 5X,'phi3',5X, 'phi4', 5X,'phi5')
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

         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K)
         DO K = 1, MS
            DO J = 1, NY
               DO I = 1, NX
                  PHI(I,J,K) = SIGMA(K)*PHI(I,J,K)
               END DO
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
          CONV = MAXVAL(ABS(PHI-PHI_OLD))/(STP*2.0D0*DT)
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
        ! In tmp_solution_file.dat, first and second columns have x & y variables,                  !
        ! respectively, density of first component followed by its phase in third and               !
        ! fourth columns, respectively, then the same for all components respectively.              !
        !                                                                                           !
        ! In convergence.dat, time, wavefunction convergence and energy convergence                 ! 
        ! parameter are written in  first and second columns, respectively.                         !
        !-------------------------------------------------------------------------------------------!

        WRITE(3,'(F14.6,6F18.10)') 2.0D0*II*DT, EN/2.0D0, (RMS(K), K = 1, MS)
        WRITE(4,'(F14.6,7F18.10)') 2.0D0*II*DT, (NORM(K), K = 1, MS), SUM(NORM), MZ
     END IF
     EN_OLD = EN
     PHI_OLD(:,:,:) = PHI(:,:,:)
     
     IF(MOD(II,NSTP).EQ.0)THEN
       WRITE (2, 1005) SUM(NORM),EN/2.0D0, (MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,NY2,K)), K = 1, MS)
       OPEN (7, FILE='./OUTPUT/tmp_solution_file_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
                STATUS='replace',FORM='FORMATTED',ACTION='WRITE')
       DO I = 1, NX
          DO J = 1, NY
             WRITE(7,'(2F8.2,10F10.5)') X(I), Y(J), (ABS(PHI(I,J,K))**2, ATAN2(AIMAG(PHI(I,J,K)),REAL(PHI(I,J,K))), K = 1, MS)
          END DO
          WRITE(7,'(11F10.2)')
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

  WRITE (2, 1006) II-1, SUM(NORM), EN/2.0D0,(MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,NY2,K)), K = 1, MS)
  1005 FORMAT('After NSTP iter.    :', 1F7.4, 6F10.5, 5F10.5)
  1006 FORMAT('After ', I8 ,' iter.:', 1F7.4, 6F10.5, 5F10.5)

  IF (SPIN.EQ.1) THEN
     WRITE (2, 1001)
  ELSE
     WRITE (2, 1011)
  END IF
  WRITE(2,*)

  !-------------------------------------------------------------------------------------------------!
  ! In solution_file.dat, first and second columns have x & y variables,                            !
  ! respectively, density of first component followed by its phase in third                         !
  ! and fourth columns, respectively, then the same for all components.                             !
  !-------------------------------------------------------------------------------------------------!

  DO I = 1, NX
     DO J = 1, NY
        WRITE(8,'(2F8.2,10F10.5)') X(I), Y(J), (ABS(PHI(I,J,K))**2, ATAN2(AIMAG(PHI(I,J,K)),REAL(PHI(I,J,K))), K = 1, MS)
     END DO
     WRITE(8,'(2F8.2,10F10.5)')
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
  1009 FORMAT("Time elapsed = ",F18.8, ' in seconds')

  CLOSE(100)
END PROGRAM CGPE2D

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Dynamic Memory Allocation to the variables and arrays                          !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE ALLOCATE_MEM()
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : NX, NY
  USE CGPE_DATA, ONLY : SPIN, X, X2, Y,Y2, V, R2, KX, KY, PHI, PHIF, PHI_OLD, SX, SY, SZ
  USE FFTW_DATA, ONLY : FFTFXY, FFTBXY

  IMPLICIT NONE
  INTEGER :: MS, STAT

  MS = 2*SPIN+1
  ALLOCATE(X(1:NX),X2(1:NX),Y(1:NY), Y2(1:NY),V(1:NX,1:NY), R2(1:NX,1:NY),KX(1:NX),KY(1:NY),STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(PHI(1:NX,1:NY,1:MS),PHIF(1:NX,1:NY,1:MS), PHI_OLD(1:NX,1:NY,1:MS),STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(SX(1:MS,1:MS),SY(1:MS,1:MS),SZ(1:MS))
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(FFTFXY(NX,NY), FFTBXY(NX,NY))
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
END SUBROUTINE ALLOCATE_MEM

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                       Memory Deallocation of the variables and arrays                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DEALLOCATE_MEM()
  USE DOUBLE_PRECISION
  USE CGPE_DATA, ONLY : X, X2, Y, Y2, V, R2, KX, KY, PHI, PHIF, PHI_OLD, SX, SY, SZ
  USE FFTW_DATA, ONLY : FFTFXY, FFTBXY

  IMPLICIT NONE
  INTEGER :: STAT

  DEALLOCATE(X, X2, Y, Y2, V, R2, KX, KY, STAT = STAT)
  DEALLOCATE(PHI, PHIF, PHI_OLD,STAT = STAT)
  DEALLOCATE(SX,SY,SZ)
  DEALLOCATE(FFTFXY, FFTBXY)
END SUBROUTINE DEALLOCATE_MEM

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Defines the SX, SY, SZ, spatial and Fourier grids, trapping potential,              !
!              and initializes the component wave functions                                        !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INITIALIZE()
  USE BASIC_DATA, ONLY : CI, PI, NX, LX, NY, LY
  USE CGPE_DATA
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY):: TMP2D, TFDENSITY
  INTEGER :: I, J
  
  SELECT CASE(SPIN)
    CASE(1)
      SX(:,:) = (SQRT(0.5D0))*RESHAPE([0.0D0,1.0D0,0.0D0,&
                1.0D0,0.0D0,1.0D0,&
                0.0D0,1.0D0,0.0D0],[3,3])
      SY(:,:) = (CI/SQRT(2.0D0))*RESHAPE([0.0D0,-1.0D0,0.0D0,&
                1.0D0,0.0D0,-1.0D0,&
                0.0D0,1.0D0,0.0D0],[3,3])
      SZ(:) = [1,0,-1]
    CASE(2)
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
  DO I = 1, 1+NX/2
     KX(I) = (I-1.0d0)*2.0D0*PI/LX
  END DO

  DO I = 1,NX/2 -1
     KX(I+1+NX/2) = -KX(1-I+NX/2)
  END DO

  !$OMP SECTION
  DO I = 1, 1+NY/2
     KY(I) = (I-1.0d0)*2.0D0*PI/LY
  END DO

  DO I = 1,NY/2 -1
     KY(I+1+NY/2)=-KY(1-I+NY/2)
  END DO

  !$OMP SECTION
  DO I = 1, NX
     X(I) = (-1.0d0 +2.0D0*(I-1.0D0)/NX)*LX/2.0D0
     X2(I) = X(I)*X(I)
  END DO

  !$OMP SECTION
  DO I = 1, NY
     Y(I) = (-1.0d0 +2.0D0*(I-1.0D0)/NY)*LY/2.0D0
     Y2(I) = Y(I)*Y(I)
  END DO
  !$OMP END PARALLEL SECTIONS 

  !$OMP PARALLEL DO PRIVATE(I,J)
  DO J = 1, NY
     DO I = 1, NX
        R2(I,J) = X2(I) + Y2(J)
        V(I,J) = ALPHAX*ALPHAX*X2(I) + ALPHAY*ALPHAY*Y2(J)
        TMP2D(I,J) = ALPHAX*X2(I)+ ALPHAY*Y2(J)
        SELECT CASE (SPIN)
          CASE(1)
            SELECT CASE (OPTION_FPC)
              CASE(1)!Initial guess wavefunctions for Ferromagnetic interactions
                IF(((4.0D0*(TAU(0)+TAU(1))/(2.0D0*PI))**(8.0D0)- ALPHAX*ALPHAX*X2(I)& 
                    - ALPHAY*ALPHAY*Y2(J)).GE.0.0D0)THEN
                  TFDENSITY(I,J) =  ((4.0D0*(TAU(0)+TAU(1))/(2.0D0*PI))**(0.5D0)&
                                    - ALPHAX*ALPHAX*X2(I) - ALPHAY*ALPHAY*Y2(J))/((TAU(0)+TAU(1)))
                ELSE
                  TFDENSITY(I,J) = 0.0D0
                END IF

                IF(((4.0D0*(TAU(0)+TAU(1))/(2.0D0*PI))**(0.25D0)).GE.10.0D0)THEN
                  PHI(I,J,1) = ((1.0D0+MAG)/2.0D0)*SQRT(TFDENSITY(I,J))
                  PHI(I,J,2) = SQRT((1.0D0-MAG**2)/2.0D0)*SQRT(TFDENSITY(I,J))
                  PHI(I,J,3) = ((1.0D0-MAG)/2.0D0)*SQRT(TFDENSITY(I,J))
                ELSE
                  PHI(I,J,1) = ((1.0D0+MAG)/2.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,2) = SQRT((1.0D0-MAG**2)/2.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,3) = ((1.0D0-MAG)/2.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                END IF

              CASE(2)!Intitial guess wavefunctions for antiferromagnetic interactions 
                PHI(I,J,1) = (SQRT((1.0D0+MAG)/2.0D0))*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                PHI(I,J,2) = 0.0D0
                PHI(I,J,3) = (SQRT((1.0D0-MAG)/2.0D0))*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))

              CASE DEFAULT!For choosing Gaussian intial guess wavefunction
                PHI(I,J,1) = EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                PHI(I,J,2) = EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                PHI(I,J,3) = EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
              END SELECT
            CASE(2)
              SELECT CASE (OPTION_FPC)
                CASE(1)!Initial guess wavefunctions for Ferromagnetic interactions
                  IF(((4.0D0*(TAU(0)+4.0D0*TAU(1))/(2.0D0*PI))**(0.5D0)&
                       - ALPHAX*ALPHAX*X2(I) - ALPHAY*ALPHAY*Y2(J)).GE.0.0D0)THEN
                    TFDENSITY(I,J) =  ((4.0D0*(TAU(0)+4.0D0*TAU(1))/(2.0D0*PI))**(0.5D0)-&
                                    ALPHAX*ALPHAX*X2(I) - ALPHAY*ALPHAY*Y2(J))/((TAU(0)+4.0D0*TAU(1)))
                  ELSE
                    TFDENSITY(I,J) = 0.0D0
                  END IF

                  IF(((4.0D0*(TAU(0)+4.0D0*TAU(1))/(2.0D0*PI))**(0.25D0)).GE.10.0D0)THEN
                    PHI(I,J,1) = (((2.0D0+MAG)**2)/16.0D0)*SQRT(TFDENSITY(I,J))
                    PHI(I,J,2) = ((SQRT(4.0D0-MAG**2)*(2.0D0+MAG))/8.0D0)*SQRT(TFDENSITY(I,J))
                    PHI(I,J,3) = ((SQRT(1.5D0)*(4.0D0-MAG**2))/8.0D0)*SQRT(TFDENSITY(I,J))
                    PHI(I,J,4) = ((SQRT(4.0D0-MAG**2)*(2.0D0-MAG))/8.0D0)*SQRT(TFDENSITY(I,J))
                    PHI(I,J,5) = (((2.0D0-MAG)**2)/16.0D0)*SQRT(TFDENSITY(I,J))
                  ELSE
                    PHI(I,J,1) = (((2.0D0+MAG)**2)/16.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                    PHI(I,J,2) = ((SQRT(4.0D0-MAG**2)*(2.0D0+MAG))/8.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                    PHI(I,J,3) = ((SQRT(1.5D0)*(4.0D0-MAG**2))/8.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                    PHI(I,J,4) = ((SQRT(4.0D0-MAG**2)*(2.0D0-MAG))/8.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                    PHI(I,J,5) = (((2.0D0-MAG)**2)/16.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                END IF

                CASE(2)!Intitial guess wavefunctions for polar interactions 
                  PHI(I,J,1) = (SQRT(2.0D0+MAG)/2.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,2) = 0.0D0
                  PHI(I,J,3) = 0.0D0
                  PHI(I,J,4) = 0.0D0
                  PHI(I,J,5) = (SQRT(2.0D0-MAG)/2.0D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))

                CASE(3)!Intitial guess wavefunctions for cyclic interactions 
                  PHI(I,J,1) = (SQRT((1.0D0+MAG)/3.0D0))*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,2) = 0.0D0
                  PHI(I,J,3) = 0.0D0
                  PHI(I,J,4) = (SQRT((2.0D0-MAG)/3.0D0))*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,5) = 0.0D0

                CASE DEFAULT !For choosing Gaussian intial guess wavefunction
                  PHI(I,J,1) = SQRT(0.2D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,2) = SQRT(0.2D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,3) = SQRT(0.2D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,4) = SQRT(0.2D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
                  PHI(I,J,5) = SQRT(0.2D0)*EXP(-TMP2D(I,J)/2.0D0)/(SQRT(PI))
              END SELECT
        END SELECT
     END DO
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
  PLANFXY = FFTW_PLAN_DFT_2D(NY,NX,FFTFXY,FFTBXY,FFTW_FORWARD,FFTW_ESTIMATE)
  PLANBXY = FFTW_PLAN_DFT_2D(NY,NX,FFTBXY,FFTFXY,FFTW_BACKWARD,FFTW_ESTIMATE)
END SUBROUTINE CREATE_PLANS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                    Destroys FFTW plans for forward and backward transforms                       !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DESTROY_PLANS()
  USE FFTW3
  USE FFTW_DATA

  IMPLICIT NONE

  CALL FFTW_DESTROY_PLAN(PLANFXY)
  CALL FFTW_DESTROY_PLAN(PLANBXY)
  CALL FFTW_CLEANUP_THREADS()
END SUBROUTINE DESTROY_PLANS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         Solves the Hamiltonian corresponding to Kinetic Energy terms in Fourier space            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KE()
  USE BASIC_DATA, ONLY : CI, NX, NY, CDT
  USE CGPE_DATA, ONLY :  PHIF, KX, KY, SPIN

  IMPLICIT NONE
  INTEGER :: I, J, K, MS

  MS = 2*SPIN+1
  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K)
  DO K = 1, MS
     DO J = 1, NY
        DO I = 1, NX
           PHIF(I,J,K) = exp(-CI*CDT*(KX(I)*KX(I) + KY(J)*KY(J)))*PHIF(I,J,K)
        END DO
     END DO
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE KE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    Solves the Hamiltonian corresponding to trapping potential and diagonal interaction terms     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SP()
  USE BASIC_DATA, ONLY : CI, NX, NY, CDT
  USE CGPE_DATA, ONLY :  V, TAU, PHI, SPIN
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  REAL(KIND=DBL), DIMENSION(1:2*SPIN+1) :: TMP
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY) :: RHO, FZ
  INTEGER :: I, J, K, MS

  MS = 2*SPIN+1

  SELECT CASE (SPIN)
    CASE(1)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,TMP,K)
      DO J = 1, NY
         DO I = 1, NX
            RHO(I,J) = SUM(ABS(PHI(I,J,1:MS))**2)
            TMP(1) = V(I,J) + TAU(0)*RHO(I,J) + TAU(1)*(ABS(PHI(I,J,1))**2 + ABS(PHI(I,J,2))**2&
                    -ABS(PHI(I,J,3))**2)
            TMP(2) = V(I,J) + TAU(0)*RHO(I,J) + TAU(1)*(ABS(PHI(I,J,1))**2 + ABS(PHI(I,J,3))**2)
            TMP(3) = V(I,J) + TAU(0)*RHO(I,J) + TAU(1)*(ABS(PHI(I,J,3))**2 + ABS(PHI(I,J,2))**2& 
                     -ABS(PHI(I,J,1))**2)
            DO K = 1, MS
               PHI(I,J,K) = PHI(I,J,K)*EXP((-CI*CDT)*TMP(K))
           END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,TMP,K)
      DO J = 1, NY
         DO I = 1, NX
            RHO(I,J) = SUM(ABS(PHI(I,J,1:MS))**2)
            FZ(I,J)  = 2.0D0*ABS(PHI(I,J,1))**2 + ABS(PHI(I,J,2))**2&
                     - ABS(PHI(I,J,4))**2 - 2.0D0*ABS(PHI(I,J,5))**2
            TMP(1)   = V(I,J) + TAU(0)*RHO(I,J) + 2.0D0*TAU(1)*FZ(I,J)&
                     + 0.4D0*TAU(2)*ABS(PHI(I,J,5))**2
            TMP(2)   = V(I,J) + TAU(0)*RHO(I,J) + TAU(1)*FZ(I,J)&
                     + 0.4D0*TAU(2)*ABS(PHI(I,J,4))**2
            TMP(3)   = V(I,J) + TAU(0)*RHO(I,J)+ 0.2D0*TAU(2)*ABS(PHI(I,J,3))**2
            TMP(4)   = (V(I,J) + TAU(0)*RHO(I,J) - TAU(1)*FZ(I,J) + 0.4D0*TAU(2)*ABS(PHI(I,J,2))**2)
            TMP(5) = (V(I,J)+TAU(0)*RHO(I,J) - 2.0D0*TAU(1)*FZ(I,J)&
                      + 0.4D0*TAU(2)*ABS(PHI(I,J,1))**2)
            DO K = 1, MS
               PHI(I,J,K) = PHI(I,J,K)*EXP((-CI*CDT)*TMP(K))
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
  END SELECT
END SUBROUTINE SP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                 Calculates FX, FY and FZ                                         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE FXYZ(PHI, FX, FY, FZ)
  USE BASIC_DATA, ONLY : NX, NY
  USE CGPE_DATA, ONLY : SPIN, SX, SY, SZ
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(OUT) :: FX, FY, FZ
  INTEGER :: I, J, MS

  MS = 2*SPIN+1

  !$OMP PARALLEL DO PRIVATE(I,J)
  DO J = 1, NY
     DO I = 1, NX
        FX(I,J) = REAL(DOT_PRODUCT(PHI(I,J,1:MS), MATMUL(SX(1:MS,1:MS),PHI(I,J,1:MS))))
        FY(I,J) = REAL(DOT_PRODUCT(PHI(I,J,1:MS), MATMUL(SY(1:MS,1:MS),PHI(I,J,1:MS))))
        FZ(I,J) = DOT_PRODUCT(SZ(1:MS),ABS(PHI(I,J,1:MS))**2)
     END DO
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE FXYZ

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            Calculates C12, C13, C23, C34, C35 and C45, i.e. the elements of Hamiltonian          !
!            corresponding to off-diagonal interaction terms in spin-2 systems                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
  USE BASIC_DATA, ONLY : NX, NY
  USE CGPE_DATA, ONLY : TAU
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:5), INTENT(IN) :: PHI
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(IN) :: FMINUS
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
  
  !$OMP PARALLEL WORKSHARE
  C12(:,:) = TAU(1)*FMINUS(:,:) - 0.4D0*TAU(2)*PHI(:,:,4)*CONJG(PHI(:,:,5)) 

  C13(:,:) = 0.2D0*TAU(2)*PHI(:,:,3)*CONJG(PHI(:,:,5))

  C23(:,:) = (SQRT(1.5D0))*TAU(1)*FMINUS(:,:)-0.2D0*TAU(2)*PHI(:,:,3)*CONJG(PHI(:,:,4)) 

  C34(:,:) = (SQRT(1.5D0))*TAU(1)*FMINUS(:,:)-0.2D0*TAU(2)*PHI(:,:,2)*CONJG(PHI(:,:,3)) 

  C35(:,:) = 0.2D0*TAU(2)*PHI(:,:,1)*CONJG(PHI(:,:,3))

  C45(:,:) = TAU(1)*FMINUS(:,:) - 0.4D0*TAU(2)*PHI(:,:,1)*CONJG(PHI(:,:,2)) 
  !$OMP END PARALLEL WORKSHARE
END SUBROUTINE MAT_C

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          Calculates the norm of individual components                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NORMC(PHI, NORM)
  USE BASIC_DATA, ONLY : NX, NY, DX, DY
  USE CGPE_DATA, ONLY : SPIN
  USE DOUBLE_PRECISION 
  
  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
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
    
  REAL(KIND=DBL), DIMENSION(1:NX) :: P
  REAL(KIND=DBL), DIMENSION(1:NY,1:2*SPIN+1) :: TMP1D
  INTEGER :: I, J, K, MS

  MS = 2*SPIN+1

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,P)
  DO K = 1, MS
     DO J = 1, NY
        DO I = 1, NX
           P(I) = ABS(PHI(I,J,K))**2
        END DO
        TMP1D(J,K) = SIMPSON(P(:), DX)
     END DO
  END DO
  !$OMP END PARALLEL DO

  DO K = 1, MS
     NORM(K) = SIMPSON(TMP1D(:,K), DY)
  END DO
END SUBROUTINE NORMC

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Calculates the total norm and normalizes the total density to 1                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NORMT(PHI, TNORM)
  USE BASIC_DATA, ONLY : NX, DX, NY, DY
  USE CGPE_DATA, ONLY : SPIN
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(INOUT) :: PHI
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

  INTEGER :: I, J, MS
  REAL(KIND=DBL), DIMENSION(1:NX) :: P
  REAL(KIND=DBL), DIMENSION(1:NY) :: TMP1D

  MS = 2*SPIN+1

  !$OMP PARALLEL DO PRIVATE(I,J,P)
  DO J = 1, NY
     DO I = 1, NX
        P(I) = SUM(ABS(PHI(I,J,1:MS))**2)
     END DO
     TMP1D(J) = SIMPSON(P(:), DX)
  END DO
  !$OMP END PARALLEL DO
  TNORM = SQRT(SIMPSON(TMP1D, DY))
  PHI = PHI/(TNORM)
END SUBROUTINE NORMT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Calculates the root mean square sizes of the (2*SPIN+1) components                  !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE RAD(PHI, RMS)
  USE BASIC_DATA, ONLY : NX, NY, DX, DY
  USE CGPE_DATA, ONLY : X2, Y2,  SPIN
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
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

  INTEGER :: I, J, K, MS
  REAL(KIND=DBL), DIMENSION(1:NX) :: TMP2D
  REAL(KIND=DBL), DIMENSION(1:NY,1:2*SPIN+1) :: TMP1D

  MS = 2*SPIN+1

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,TMP2D)
  DO K = 1, MS
     DO J = 1, NY
        DO I = 1, NX
           TMP2D(I) = (X2(I)+Y2(J))*ABS(PHI(I,J,K))**2
        END DO
        TMP1D(J,K) = SIMPSON(TMP2D(:), DX)
     END DO
  END DO
  !$OMP END PARALLEL DO

  DO J = 1, MS
     RMS(J) = SQRT(SIMPSON(TMP1D(:,J), DY))
  END DO
END SUBROUTINE RAD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       Calculates the (2*SPIN+1) component chemical potentials, energy, and magnetization         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE ENERGY(PHI, MU, EN, MZ)
  USE BASIC_DATA, ONLY : NX, NY, DX, DY, CI, EPS
  USE CGPE_DATA, ONLY : V, TAU, SPIN
  USE SOC_DATA, ONLY : GAMMAX, GAMMAY
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
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
      USE BASIC_DATA, ONLY : NX, NY, CI
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(OUT) :: FX, FY, FZ
    END SUBROUTINE FXYZ

    SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      USE BASIC_DATA, ONLY : NX, NY
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:5), INTENT(IN) :: PHI
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(IN) ::  FMINUS
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
    END SUBROUTINE MAT_C

    SUBROUTINE NORMC(PHI, NORM)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: NORM
    END SUBROUTINE NORMC
  END INTERFACE
!---------------------------------------------------------------------------------------------------

  INTEGER :: I, J, K
  REAL(KIND=DBL), DIMENSION(1:NX,2*SPIN+3) :: TMP2D
  REAL(KIND=DBL), DIMENSION(1:NX, 1:NY) :: FX, FY, FZ, RHO
  REAL(KIND=DBL), DIMENSION(1:NY,1:2*SPIN+3) :: TMP1D
  COMPLEX(KIND=DBL), DIMENSION(1:NX, 1:NY,1:2*SPIN+1) :: DPX, DPY
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:2*SPIN+1) :: DP2, DPX2
  COMPLEX(KIND=DBL), DIMENSION(1:NX, 1:NY) :: THETA, FMINUS, C12, C13, C23, C34, C35, C45
  REAL(KIND=DBL), DIMENSION(1:2*SPIN+1) :: NORM
 
  CALL FXYZ(PHI, FX, FY, FZ)
  
  CALL NORMC(PHI, NORM)

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(K,J)
  DO K = 1, 2*SPIN+1
     DO J = 1, NY
        DPX2(:,J,K) = DIFF(REAL(PHI(:,J,K)), DX)**2 + DIFF(AIMAG(PHI(:,J,K)), DX)**2
     END DO
  END DO
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(K,I)
  DO K = 1, 2*SPIN+1
     DO I = 1, NX
        DP2(I,:,K) = DPX2(I,:,K) + DIFF(REAL(PHI(I,:,K)), DY)**2 + DIFF(AIMAG(PHI(I,:,K)), DY)**2
        DPY(I,:,K) = DIFF(REAL(PHI(I,:,K)), DY) + CI*DIFF(AIMAG(PHI(I,:,K)), DY)
     END DO
  END DO
 !$OMP END PARALLEL DO

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(K,J) 
  DO K = 1, 2*SPIN+1
     DO J = 1, NY
        DPX(:,J,K) = DIFF(REAL(PHI(:,J,K)), DX) + CI*DIFF(AIMAG(PHI(:,J,K)), DX)
     END DO
  END DO
  !$OMP END PARALLEL DO

  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO PRIVATE(I,J,TMP2D)
      DO J = 1, NY
         DO I = 1, NX
            RHO(I,J) = SUM(ABS(PHI(I,J,1:2*SPIN+1))**2)
            TMP2D(I,1) = REAL(V(I,J)*RHO(I,J) + TAU(0)*RHO(I,J)**2/2.0D0 +&
                         TAU(1)*(ABS(FX(I,J))**2 + ABS(FY(I,J))**2 + ABS(FZ(I,J))**2)/2.0D0+&
                         DP2(I,J,1) + DP2(I,J,2) + DP2(I,J,3) +&
                         SQRT(2.0D0)*(CONJG(PHI(I,J,1))*(-CI*GAMMAX*DPX(I,J,2)&
                         -GAMMAY*DPY(I,J,2))+CONJG(PHI(I,J,2))*(-CI*GAMMAX*DPX(I,J,3)-GAMMAY*DPY(I,J,3))+&
                         CONJG(PHI(I,J,2))*(-CI*GAMMAX*DPX(I,J,1) + GAMMAY*DPY(I,J,1)) &
                         +CONJG(PHI(I,J,3))*(-CI*GAMMAX*DPX(I,J,2) + GAMMAY*DPY(I,J,2))))

            TMP2D(I,2) = ABS(PHI(I,J,1))**2-ABS(PHI(I,J,3))**2

            TMP2D(I,3) = REAL(V(I,J)*ABS(PHI(I,J,1))*ABS(PHI(I,J,1)) + DP2(I,J,1) +&
                         TAU(0)*(ABS(PHI(I,J,1))**2+ABS(PHI(I,J,2))**2+ABS(PHI(I,J,3))**2)*ABS(PHI(I,J,1))**2+&
                         TAU(1)*(ABS(PHI(I,J,1))**2+ABS(PHI(I,J,2))**2-ABS(PHI(I,J,3))**2)*ABS(PHI(I,J,1))**2+&
                         TAU(1)*CONJG(PHI(I,J,3))*PHI(I,J,2)**2*CONJG(PHI(I,J,1))+&
                         SQRT(2.0D0)*(CONJG(PHI(I,J,1))*(-CI*GAMMAX*DPX(I,J,2)-GAMMAY*DPY(I,J,2))))

            TMP2D(I,4) = REAL(V(I,J)*ABS(PHI(I,J,2))*ABS(PHI(I,J,2)) + DP2(I,J,2) +&
                         TAU(0)*(ABS(PHI(I,J,1))**2+ABS(PHI(I,J,2))**2+ABS(PHI(I,J,3))**2)*ABS(PHI(I,J,2))**2+&
                         TAU(1)*(ABS(PHI(I,J,1))**2+ABS(PHI(I,J,3))**2)*ABS(PHI(I,J,2))**2+&
                         2.0D0*TAU(1)*PHI(I,J,1)*PHI(I,J,3)*CONJG(PHI(I,J,2))**2+&
                         SQRT(2.0D0)*(CONJG(PHI(I,J,2))*(-CI*GAMMAX*DPX(I,J,3) - GAMMAY*DPY(I,J,3))+&
                         CONJG(PHI(I,J,2))*(-CI*GAMMAX*DPX(I,J,1) + GAMMAY*DPY(I,J,1))))

            TMP2D(I,5) = REAL(V(I,J)*ABS(PHI(I,J,3))*ABS(PHI(I,J,3)) + DP2(I,J,3) +&
                         TAU(0)*(ABS(PHI(I,J,1))**2+ABS(PHI(I,J,2))**2+ABS(PHI(I,J,3))**2)*ABS(PHI(I,J,3))**2+&
                         TAU(1)*(ABS(PHI(I,J,2))**2+ABS(PHI(I,J,3))**2-ABS(PHI(I,J,1))**2)*ABS(PHI(I,J,3))**2+&
                         TAU(1)*CONJG(PHI(I,J,1))*PHI(I,J,2)**2*CONJG(PHI(I,J,3))+&
                         SQRT(2.0D0)*(CONJG(PHI(I,J,3))*(-CI*GAMMAX*DPX(I,J,2)+GAMMAY*DPY(I,J,2))))

         END DO
         DO I = 1, 2*SPIN+3
            TMP1D(J,I) = SIMPSON(TMP2D(:,I), DX)
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      FMINUS = FX - CI * FY
      CALL MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      !$OMP PARALLEL DO PRIVATE(I,J,TMP2D)
      DO J = 1, NY
         DO I = 1, NX
            RHO(I,J) = SUM(ABS(PHI(I,J,1:2*SPIN+1))**2)
            THETA(I,J) = (2.0D0*PHI(I,J,1)*PHI(I,J,5)-2.0D0*PHI(I,J,2)*PHI(I,J,4)+PHI(I,J,3)**2)/SQRT(5.0D0)
            TMP2D(I,1) = REAL(V(I,J)*RHO(I,J) + TAU(0)*RHO(I,J)**2/2.0D0 +&
                         TAU(1)*(ABS(FX(I,J))**2 + ABS(FY(I,J))**2 + ABS(FZ(I,J))**2)/2.0D0+&
                         TAU(2)*(THETA(I,J)*CONJG(THETA(I,J)))/2.0D0 + DP2(I,J,1) + DP2(I,J,2)+&
                         DP2(I,J,3)+ DP2(I,J,4)+ DP2(I,J,5)&
                         -CI*2.0D0*(CONJG(PHI(I,J,1))*GAMMAX*DPX(I,J,2)-CI*CONJG(PHI(I,J,1))*GAMMAY*DPY(I,J,2) +&
                         CONJG(PHI(I,J,2))*GAMMAX*DPX(I,J,1)+CI*CONJG(PHI(I,J,2))*GAMMAY*DPY(I,J,1)+SQRT(3.0D0/2.0D0)*&
                         (CONJG(PHI(I,J,2))*GAMMAX*DPX(I,J,3) -CI*CONJG(PHI(I,J,2))*GAMMAY*DPY(I,J,3) +&
                         CONJG(PHI(I,J,3))*GAMMAX*DPX(I,J,2)+&
                         CI*CONJG(PHI(I,J,3))*GAMMAY*DPY(I,J,2)+CONJG(PHI(I,J,3))*GAMMAX*DPX(I,J,4)-&
                         CI*CONJG(PHI(I,J,3))*GAMMAY*DPY(I,J,4)+&
                         CONJG(PHI(I,J,4))*GAMMAX*DPX(I,J,3)+CI*CONJG(PHI(I,J,4))*GAMMAY*DPY(I,J,3))+&
                         CONJG(PHI(I,J,4))*GAMMAX*DPX(I,J,5)-&
                         CI*CONJG(PHI(I,J,4))*GAMMAY*DPY(I,J,5)+CONJG(PHI(I,J,5))*GAMMAX*DPX(I,J,4)+&
                         CI*CONJG(PHI(I,J,5))*GAMMAY*DPY(I,J,4)))
 
            TMP2D(I,2) = 2.0D0*ABS(PHI(I,J,1))**2 + ABS(PHI(I,J,2))**2 - ABS(PHI(I,J,4))**2 - &
                         2.0D0*ABS(PHI(I,J,5))**2

            TMP2D(I,3) = REAL((V(I,J)+TAU(0)*RHO(I,J) + 2.0D0*TAU(1)*FZ(I,J)+&
                         0.4D0*TAU(2)*ABS(PHI(I,J,5))**2)*ABS(PHI(I,J,1))**2+&
                         (C12(I,J)*PHI(I,J,2)+C13(I,J)*PHI(I,J,3))*CONJG(PHI(I,J,1))&
                         -CI*2.0D0*(CONJG(PHI(I,J,1))*GAMMAX*DPX(I,J,2)-CI*CONJG(PHI(I,J,1))*GAMMAY*DPY(I,J,2))+&
                         ABS(DPX(I,J,1))**2+ ABS(DPY(I,J,1))**2)

            TMP2D(I,4) = REAL((V(I,J)+TAU(0)*RHO(I,J) + TAU(1)*FZ(I,J)+&
                         0.4D0*TAU(2)*ABS(PHI(I,J,4))**2)*ABS(PHI(I,J,2))**2+&
                         (CONJG(C12(I,J))*PHI(I,J,1)+C23(I,J)*PHI(I,J,3))*CONJG(PHI(I,J,2))&
                         -CI*2.0D0*(CONJG(PHI(I,J,2))*GAMMAX*DPX(I,J,1)+CI*CONJG(PHI(I,J,2))*GAMMAY*DPY(I,J,1)+&
                         SQRT(1.5D0)*(CONJG(PHI(I,J,2))*GAMMAX*DPX(I,J,3) -CI*CONJG(PHI(I,J,2))*GAMMAY*DPY(I,J,3)))+&
                         ABS(DPX(I,J,2))**2+ ABS(DPY(I,J,2))**2)

            TMP2D(I,5) = REAL((V(I,J)+TAU(0)*RHO(I,J) +&
                         0.2D0*TAU(2)*ABS(PHI(I,J,3))**2)*ABS(PHI(I,J,3))**2+&
                         (CONJG(C13(I,J))*PHI(I,J,1)+CONJG(C23(I,J))*PHI(I,J,2) &
                         + C34(I,J)*PHI(I,J,4)+C35(I,J)*PHI(I,J,5))*CONJG(PHI(I,J,3))&
                         -CI*2.0D0*(SQRT(1.5D0)*(CONJG(PHI(I,J,3))*GAMMAX*DPX(I,J,2) +CI*CONJG(PHI(I,J,3))*GAMMAY*DPY(I,J,2)+& 
                         CONJG(PHI(I,J,3))*GAMMAX*DPX(I,J,4)-CI*CONJG(PHI(I,J,3))*GAMMAY*DPY(I,J,4)))+&
                         ABS(DPX(I,J,3))**2+ ABS(DPY(I,J,3))**2)

            TMP2D(I,6) = REAL((V(I,J)+TAU(0)*RHO(I,J) - TAU(1)*FZ(I,J)+&
                         0.4D0*TAU(2)*ABS(PHI(I,J,2))**2)*ABS(PHI(I,J,4))**2+&
                         (CONJG(C34(I,J))*PHI(I,J,3)+C45(I,J)*PHI(I,J,5))*CONJG(PHI(I,J,4))&
                         -CI*2.0D0*(CONJG(PHI(I,J,4))*GAMMAX*DPX(I,J,5)-CI*CONJG(PHI(I,J,4))*GAMMAY*DPY(I,J,5)+&
                         SQRT(1.5D0)*(CONJG(PHI(I,J,4))*GAMMAX*DPX(I,J,3) + CI*CONJG(PHI(I,J,4))*GAMMAY*DPY(I,J,3)))+&
                         ABS(DPX(I,J,4))**2+ ABS(DPY(I,J,4))**2)
   
            TMP2D(I,7) = REAL((V(I,J)+TAU(0)*RHO(I,J) - 2.0D0*TAU(1)*FZ(I,J)+&
                         0.4D0*TAU(2)*ABS(PHI(I,J,1))**2)*ABS(PHI(I,J,5))**2+&
                         (CONJG(C35(I,J))*PHI(I,J,3)+CONJG(C45(I,J))*PHI(I,J,4))*CONJG(PHI(I,J,5))&
                         -CI*2.0D0*(CONJG(PHI(I,J,5))*GAMMAX*DPX(I,J,4)+CI*CONJG(PHI(I,J,5))*GAMMAY*DPY(I,J,4))+&
                         ABS(DPX(I,J,5))**2+ ABS(DPY(I,J,5))**2)
         END DO
         DO I = 1, 2*SPIN+3
            TMP1D(J,I) = SIMPSON(TMP2D(:,I), DX)
         END DO
      END DO
      !$OMP END PARALLEL DO 
  END SELECT

  EN = SIMPSON(TMP1D(:,1),DY)/SUM(NORM)
  MZ = SIMPSON(TMP1D(:,2),DY)
  DO K = 1, 2*SPIN+1
     MU(K) = SIMPSON(TMP1D(:,K+2),DY)/(NORM(K)+ EPS)
  END DO
END SUBROUTINE ENERGY 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                  Solves the non-linear equations by using Newton Raphson method                  !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NEWTON(X, N1, N2, N3, N4, N5, MAG, ERR_MSG)
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : TOL_NR
  USE CGPE_DATA, ONLY : NTRIAL

  IMPLICIT NONE

  REAL(KIND=DBL), INTENT(IN) :: N1, N2, N3, N4, N5, MAG
  INTEGER, INTENT(INOUT) :: ERR_MSG
  REAL(KIND=DBL), DIMENSION(1:2), INTENT(INOUT) :: X

  !USES usrfun
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
     FFTFXY(:,:) = PHI(:,:,L)
     CALL FFTW_EXECUTE_DFT(PLANFXY,FFTFXY,FFTBXY)
     PHIF(:,:,L) = FFTBXY(:,:)
  END DO
END SUBROUTINE FFT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         Calculates the discrete backward Fourier transform of component wavefunctions            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE BFT()
  USE BASIC_DATA, ONLY : NX, NY
  USE CGPE_DATA,  ONLY : PHI, PHIF, SPIN
  USE FFTW3, ONLY : FFTW_EXECUTE_DFT
  USE FFTW_DATA

  IMPLICIT NONE
  INTEGER :: L, MS

  MS = 2*SPIN+1
  DO L = 1, MS
     FFTBXY(:,:) = PHIF(:,:,L)
     CALL FFTW_EXECUTE_DFT(PLANBXY,FFTBXY,FFTFXY)
     PHI(:,:,L) = FFTFXY(:,:)/DBLE(NX*NY)
  END DO
END SUBROUTINE BFT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               Solves the Hamiltonian corresponding to off-diagonal interaction terms             !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SE()
  USE BASIC_DATA, ONLY : CI, NX, NY, CDT, EPS
  USE CGPE_DATA, ONLY : PHI, TAU, SPIN
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE FXYZ(PHI, FX, FY, FZ)
      USE BASIC_DATA, ONLY : NX, NY, CI
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(OUT) :: FX, FY, FZ
    END SUBROUTINE FXYZ

    SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      USE BASIC_DATA, ONLY : NX, NY
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:5), INTENT(IN) :: PHI
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(IN) :: FMINUS
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
    END SUBROUTINE MAT_C
  END INTERFACE

  INTEGER :: I, J, K
  COMPLEX(KIND=DBL) :: A1, B1, OMEGA, DELTA, KAPPA
  COMPLEX(KIND=DBL), DIMENSION(1:3) :: TMP

  REAL(KIND=DBL), DIMENSION(1:NX,1:NY) :: FX, FY, FZ, FX1, FY1, FZ1
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY) :: FMINUS,&
                                 C12, C13, C23, C34, C35, C45,&
                                 C121, C131, C231, C341, C351, C451,&
                                 FMINUS1
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:5) :: PHI1

  EXTERNAL :: ZHEEV
  INTEGER, PARAMETER :: LDA = 5
  INTEGER, PARAMETER :: LWMAX = 1000
  INTEGER :: INFO, LWORK
  REAL(KIND=DBL) :: RWORK(13), W(5)
  COMPLEX(KIND=DBL) :: A(1:5,1:5), WORK(LWMAX)

  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,TMP,A1,B1,OMEGA,DELTA,KAPPA)
      DO J = 1, NY
         DO I = 1, NX
            A1 =  (TAU(1) * PHI(I,J,2)*  CONJG(PHI(I,J,3)))
            B1 =  (TAU(1) * PHI(I,J,2)*  CONJG(PHI(I,J,1)))
            OMEGA = CDT*SQRT((ABS(A1))**2+(ABS(B1))**2)
  
            IF(ABS(OMEGA).GT.EPS)THEN
               DELTA = (COS(OMEGA) -1)/(OMEGA**2)
               KAPPA = -CI * SIN(OMEGA)/OMEGA
            ELSE
               DELTA = -0.5D0
               KAPPA = -CI
            END IF
  
            TMP(1) = (DELTA*CDT**2)*(A1*CONJG(A1)*PHI(I,J,1)+&
                     A1*CONJG(B1)*PHI(I,J,3)) + KAPPA*CDT*A1&
                    *PHI(I,J,2)
            TMP(2) = (DELTA*CDT**2)*(A1*CONJG(A1)+B1*CONJG(B1))&
                    *PHI(I,J,2) + KAPPA*CDT*&
                     (CONJG(A1)*PHI(I,J,1)+CONJG(B1)*PHI(I,J,3))
            TMP(3) = (DELTA*CDT**2)*(CONJG(A1)*B1*&
                     PHI(I,J,1)+B1*CONJG(B1)*PHI(I,J,3)) + KAPPA*CDT*&
                     B1*PHI(I,J,2)
  
            PHI(I,J,1) = PHI(I,J,1) + TMP(1)
            PHI(I,J,2) = PHI(I,J,2) + TMP(2)
            PHI(I,J,3) = PHI(I,J,3) + TMP(3)
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      CALL FXYZ(PHI, FX, FY,FZ)
   
      FMINUS(:,:) = FX(:,:) - CI * FY(:,:)
   
      CALL MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
   
      !$OMP PARALLEL WORKSHARE
      PHI1(:,:,1) = PHI(:,:,1) - (CI*CDT)*(C12(:,:)*PHI(:,:,2)+C13(:,:)*PHI(:,:,3))
      PHI1(:,:,2) = PHI(:,:,2) - (CI*CDT)*(CONJG(C12(:,:))*PHI(:,:,1)+C23(:,:)*PHI(:,:,3))
      PHI1(:,:,3) = PHI(:,:,3) - (CI*CDT)*(CONJG(C13(:,:))*PHI(:,:,1)&
                    +CONJG(C23(:,:))*PHI(:,:,2)+C34(:,:)*PHI(:,:,4)+C35(:,:)*PHI(:,:,5))
      PHI1(:,:,4) = PHI(:,:,4) - (CI*CDT)*(CONJG(C34(:,:))*PHI(:,:,3)+C45(:,:)*PHI(:,:,5))
      PHI1(:,:,5) = PHI(:,:,5) - (CI*CDT)*(CONJG(C35(:,:))*PHI(:,:,3)+CONJG(C45(:,:))*PHI(:,:,4))
      !$OMP END PARALLEL WORKSHARE
    
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
    
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,W,LWORK,RWORK,INFO,WORK) FIRSTPRIVATE(A)
      DO J = 1, NY
        DO I = 1, NX
           A(1,2) = C12(I,J)
           A(1,3) = C13(I,J)
   
           A(2,1) = CONJG(C12(I,J))
           A(2,3) = C23(I,J)
   
           A(3,1) = CONJG(C13(I,J))
           A(3,2) = CONJG(C23(I,J))
   
           A(3,4) = C34(I,J)
           A(3,5) = C35(I,J)
   
           A(4,3) = CONJG(C34(I,J))
           A(4,5) = C45(I,J)
   
           A(5,3) = CONJG(C35(I,J))
           A(5,4) = CONJG(C45(I,J))
   
           LWORK = -1
           CALL ZHEEV( 'Vectors', 'Lower', 5, A, LDA, W, WORK, LWORK, RWORK, INFO )
           LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
           CALL ZHEEV( 'Vectors', 'Lower', 5, A, LDA, W, WORK, LWORK, RWORK, INFO )
   
           PHI(I,J,1:5) = MATMUL(CONJG(TRANSPOSE(A(1:5,1:5))), PHI(I,J,1:5))
           DO K = 1, 5
              PHI(I,J,K) = EXP(-CI*CDT*W(K)) * PHI(I,J,K)
           END DO
           PHI(I,J,1:5) = MATMUL(A(1:5,1:5), PHI(I,J,1:5))
   
           A = CMPLX(0.0D0,0.0D0,KIND=DBL)
           W = 0.0D0
        END DO
      END DO
      !$OMP END PARALLEL DO
    END SELECT
END SUBROUTINE SE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       Solves the Hamiltonian corresponding to spin-orbit coupling terms in Fourier space         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SOC()
  USE BASIC_DATA, ONLY : CI, NX, NY, EPS, CDT
  USE CGPE_DATA, ONLY : KX, KY, PHIF, SPIN
  USE SOC_DATA, ONLY : GAMMAX, GAMMAY
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  INTEGER :: I, J, K, MS
  COMPLEX(KIND=DBL), DIMENSION(1:2*SPIN+1,1:2*SPIN+1) :: A
  REAL(KIND=DBL) :: ALPHA, BETA

  MS = 2*SPIN+1
  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,ALPHA,BETA,A)
      DO J = 1, NY
         DO I = 1, NX
            ALPHA= 2.0D0*SQRT((GAMMAX**2)*KX(I)**2+(GAMMAY**2)*KY(J)**2)
            IF(ALPHA.LE.EPS)THEN 
              BETA = 0.0D0 
            ELSE
              BETA = ATAN2(GAMMAY*KY(J),GAMMAX*KX(I))
            END IF   

            A(1,1) = 0.5D0*EXP(-2.0D0*CI*BETA)
            A(2,1) = -SQRT(0.5D0)*EXP(-CI*BETA)
            A(3,1) = 0.5D0

            A(1,2) = -SQRT(0.5D0)*EXP(-2.0D0*CI*BETA)
            A(2,2) = 0.0d0
            A(3,2) = SQRT(0.5D0)

            A(1,3) = A(1,1)
            A(2,3) = -A(2,1)
            A(3,3) = A(3,1)

            PHIF(I,J,1:MS) = MATMUL(CONJG(TRANSPOSE(A(1:MS,1:MS))), PHIF(I,J,1:MS))

            DO K = -SPIN, SPIN, 1
               PHIF(I,J,K+SPIN+1) = EXP(-CI*CDT*K*ALPHA)*PHIF(I,J,K+SPIN+1)
            END DO
            PHIF(I,J,1:MS) = MATMUL(A(1:MS,1:MS), PHIF(I,J,1:MS))
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,ALPHA,BETA,A)
      DO J = 1, NY
         DO I = 1, NX
            ALPHA= 2.0D0*SQRT((GAMMAX**2)*KX(I)**2+(GAMMAY**2)*KY(J)**2)
            IF(ALPHA.LE.EPS)THEN
              BETA = 0.0D0
            ELSE
              BETA = ATAN2(GAMMAY*KY(J),GAMMAX*KX(I))
            END IF

            A(1,1) = 0.25D0*EXP(-2.0D0*CI*BETA)
            A(2,1) = -0.5D0*EXP(-CI*BETA)
            A(3,1) = SQRT(0.375D0)
            A(4,1) = CONJG(A(2,1))
            A(5,1) = CONJG(A(1,1))
     
            A(1,2) = -0.5D0*EXP(-2.0D0*CI*BETA)
            A(2,2) = 0.5D0*EXP(-CI*BETA)
            A(3,2) = CMPLX(0.0D0,0.0D0,KIND=DBL)
            A(4,2) = -CONJG(A(2,2))
            A(5,2) = -CONJG(A(1,2))

            A(1,3) = SQRT(0.375D0)*EXP(-2.0D0*CI*BETA)
            A(2,3) = A(3,2)
            A(3,3) = -CMPLX(0.5D0,0.0D0,KIND=DBL)
            A(4,3) = A(3,2)
            A(5,3) = CONJG(A(1,3))
     
            A(1,4) = -0.5D0*EXP(-2.0D0*CI*BETA)
            A(2,4) = -0.5D0*EXP(-CI*BETA)
            A(3,4) = A(3,2)!0.0D0
            A(4,4) = -CONJG(A(2,4))
            A(5,4) = -CONJG(A(1,4))

            A(1,5) = 0.25D0*EXP(-2.0D0*CI*BETA)
            A(2,5) = 0.5D0*EXP(-CI*BETA)
            A(3,5) = SQRT(0.375D0)
            A(4,5) = CONJG(A(2,5))
            A(5,5) = CONJG(A(1,5))

            PHIF(I,J,1:MS) = MATMUL(CONJG(TRANSPOSE(A(1:MS,1:MS))), PHIF(I,J,1:MS))
      
            DO K = -SPIN, SPIN, 1
               PHIF(I,J,K+SPIN+1) = EXP(-CI*CDT*K*ALPHA)*PHIF(I,J,K+SPIN+1)
            END DO
            PHIF(I,J,1:MS) = MATMUL(A(1:MS,1:MS), PHIF(I,J,1:MS))     
         END DO
      END DO
      !$OMP END PARALLEL DO
  END SELECT
END SUBROUTINE SOC
 
























