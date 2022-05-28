! File name: cgpe3d.f90
! Title: FORTRESS: FORTRAN programs to solve coupled
! Gross-Pitaevskii equations for spin-orbit coupled spin-f
! Bose-Einstein condensate with spin f = 1 or 2
!        
! Authors: Paramjeet Banger, Pardeep Kaur, Arko Roy, & Sandeep Gautam
!
! Fortran code to solve 3D coupled Gross-Pitaevskii equations (CGPEs)
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
! NZ = Number of space grid points in z-direction
! DX = Spatial step size along x
! DY = Spatial step size along y
! DZ = Spatial step size along z
! DT = Temporal step size
! LX = NX*DX, LY = NY*DY, LZ = NZ*DZ (3D spatial domain)
! AMU = atomic mass unit in Kg
! HBAR = reduced Planck's constat in SI units
! CDT = Complex varible which will be defined as '-idt' and 'dt' in
!       imaginary and realtime propagations, respectively 
!---------------------------------------------------------------------------------------------------

MODULE DOUBLE_PRECISION
INTEGER, PARAMETER :: DBL=KIND(0.0D0)
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
  INTEGER :: NX, NX2, NXX
  INTEGER :: NY, NY2, NYY
  INTEGER :: NZ, NZ2, NZZ
  REAL(KIND=DBL) :: DX, DY, DZ, DT
  REAL(KIND=DBL) :: LX, LY, LZ
  INTEGER :: STP, NSTP
  COMPLEX(KIND=DBL) :: CDT
END MODULE BASIC_DATA


MODULE CGPE_DATA
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : NX, NY,NZ
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
  INTEGER :: SPIN, NATOMS
  INTEGER, PARAMETER :: NTRIAL = 1000
  REAL(KIND=DBL) :: TOL
  INTEGER :: SWITCH_IM, OPTION_FPC
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: X, X2
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: Y, Y2
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: Z, Z2
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: KX
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: KY
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: KZ
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:) :: V, R2
  COMPLEX(KIND=DBL), ALLOCATABLE, DIMENSION(:,:,:,:) :: PHI, PHIF, PHI_OLD
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
  REAL(KIND=DBL):: GAMMAX, GAMMAY, GAMMAZ
  INTEGER :: SWITCH_SOC 
END MODULE SOC_DATA

MODULE FFTW_DATA
  USE BASIC_DATA, ONLY : NX, NY, NZ
  USE FFTW3, ONLY : C_DOUBLE_COMPLEX, C_PTR, C_INT
  COMPLEX(KIND=C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:) :: FFTFXYZ, FFTBXYZ
  TYPE(C_PTR) :: PLANFX, PLANBX, PLANFXYZ, PLANBXYZ
  INTEGER(KIND=C_INT) :: THREADS_INIT
END MODULE FFTW_DATA
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                       Main Program - CGPE3D                                      ! 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PROGRAM CGPE3D
  USE BASIC_DATA
  USE CGPE_DATA
  USE OMP_LIB
  USE SOC_DATA
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
      COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: NORM
    END SUBROUTINE NORMC
 
    SUBROUTINE NORMT(PHI, TNORM)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(INOUT) :: PHI
      REAL(KIND=DBL), INTENT(OUT) :: TNORM
    END SUBROUTINE NORMT

    SUBROUTINE RAD(PHI, RMS)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: RMS
    END SUBROUTINE RAD

    SUBROUTINE ENERGY(PHI, MU, EN, MZ)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
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

  INTEGER :: I, J, K, L, II, MS, IOSTAT, ERR_MSG
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: NORM, RMS, MU, SIGMA
  REAL(KIND=DBL) :: EN, EN_OLD, MZ, START, FINISH
  REAl(KIND=DBL), DIMENSION(1:2) :: SIG
  REAL(KIND=DBL) :: CONV = 0.0D0, CONV_EN=0.0D0
  REAL(KIND=DBL) :: TNORM
  REAL(KIND=DBL), ALLOCATABLE, DIMENSION(:) :: TMP
  CHARACTER (LEN=10) DATE, TIME, ZONE
  INTEGER, DIMENSION(1:8) :: VALUES
  CHARACTER(LEN=1) :: SPIN_ID
  LOGICAL :: FEXIST

  !$ START = OMP_GET_WTIME()
  CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES)
  
  OPEN (100, FILE='err_msg.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
              ACTION='WRITE')

  INQUIRE (FILE='./INPUT/input3D.dat', EXIST = FEXIST)
  IF(FEXIST)THEN
     OPEN (10, FILE='./INPUT/input3D.dat', STATUS='OLD',FORM='FORMATTED',&
              ACTION='READ')
  ELSE
    WRITE(100,*) "'input3D.dat' does not exist."
    PRINT*, "'input3D.dat' does not exist."
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
  READ(10,*) NX, NY, NZ
  READ(10,*) DX, DY, DZ, DT
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
  READ(10,*) GAMMAX, GAMMAY, GAMMAZ
  IF(ABS(GAMMAX).GT.EPS .OR. ABS(GAMMAY).GT.EPS .OR. ABS(GAMMAZ).GT.EPS) SWITCH_SOC = 1
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
  NZ2 = NZ/2+1
  NZZ = NZ-1

  ALPHAX = NUX/NUX
  ALPHAY = NUY/NUX
  ALPHAZ = NUZ/NUX

  LX = NX*DX
  LY = NY*DY
  LZ = NZ*DZ

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
      OPEN (2, FILE='./OUTPUT/file1_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
      OPEN (3, FILE='./OUTPUT/file2_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
      OPEN (4, FILE='./OUTPUT/file3_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
      OPEN (8, FILE='./OUTPUT/solution_file_im_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
      OPEN (9, FILE='./OUTPUT/convergence.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
    CASE(0)
      CDT = CMPLX(DT,0.0D0,KIND=DBL)
      OPEN (1, FILE='./INPUT/initial_sol_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', IOSTAT=iostat, STATUS='OLD',FORM='FORMATTED',&
           ACTION='READ')
      OPEN (2, FILE='./OUTPUT/file1_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
      OPEN (3, FILE='./OUTPUT/file2_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
      OPEN (4, FILE='./OUTPUT/file3_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
      OPEN (8, FILE='./OUTPUT/solution_file_re_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat', STATUS='UNKNOWN',FORM='FORMATTED',&
           ACTION='WRITE')
  END SELECT

  WRITE(2,1007) VALUES(3),VALUES(2),VALUES(1),VALUES(5),VALUES(6),VALUES(7)

  OMEGAM = 2.0D0*PI*NUX
  AOSC = SQRT(HBAR/(M*AMU*OMEGAM))
  
  SELECT CASE(SPIN)
    CASE(1)
      TAU(0) = 2.0D0*(4.0D0*PI*REAL(NATOMS)*(A0*AU+2.0D0*A2*AU)/(3.0D0*AOSC))
      TAU(1) = 2.0D0*(4.0D0*PI*REAL(NATOMS)*(A2*AU-A0*AU)/(3.0D0*AOSC))
    CASE(2)
      TAU(0) = 2.0D0*(4.0D0*PI*DBLE(NATOMS)*(4.0D0*A2*AU+3.0D0*A4*AU)/(7.0D0*AOSC))
      TAU(1) = 2.0D0*(4.0D0*PI*DBLE(NATOMS)*(A4*AU-A2*AU)/(7.0D0*AOSC))
      TAU(2) = 2.0D0*(4.0D0*PI*DBLE(NATOMS)*(7.0D0*A0*AU-10.0D0*A2*AU+3.0D0*A4*AU)/(7.0D0*AOSC))
  END SELECT
  WRITE(2,*)

  WRITE(2,899) OPENMP_THREADS, FFTW_THREADS
  WRITE(2,900) SWITCH_IM, OPTION_FPC, SWITCH_SOC, GAMMAX, GAMMAY, GAMMAZ, TOL
  WRITE(2,901) ALPHAX, ALPHAY, ALPHAZ
  WRITE(2,*)
  WRITE(2,902) NX, NY, NZ
  WRITE(2,903) NITER, NSTP, STP
  WRITE(2,904) DX, DY, DZ
  WRITE(2,905) DT, MAG
  WRITE(2,*)
  WRITE(2,908) AOSC
  IF (SPIN .EQ. 1) THEN
     WRITE(2,906) A0, A2
     WRITE(2,907) TAU(0)/2.0D0, TAU(1)/2.0D0
  ELSE
     WRITE(2,911) A0, A2, A4
     WRITE(2,912) TAU(0)/2.0D0, TAU(1)/2.0D0, TAU(2)/2.0D0
  END IF
  WRITE(2,*)

  899 FORMAT(' OPENMP_THREADS = ', I3,', FFTW_THREADS = ', I3)
  900 FORMAT(' SWITCH_IM = ', I2,', OPTION_FPC = ', I2, ', SWITCH_SOC = ',&
                I2, ', GAMMAX = ', F6.2, ', GAMMAY = ', F6.2, ', GAMMAZ = ', F6.2, ', TOL = ', ES9.2)
  901 FORMAT(' Anisotropy ALPHAX =',F12.6,', ALPHAY = ',F12.6, ', ALPHAZ = ',F12.6)
  902 FORMAT(' No. of space steps NX = ',I8, ', NY = ', I8, ', NZ = ', I8)
  903 FORMAT(' No. of time steps : NITER = ',I9,', NSTP = ',I9,', STP = ',I9)
  904 FORMAT(' Spatial step-size DX = ',F10.6,', DY = ',F10.6,', DZ = ',F10.6)
  905 FORMAT(' Temporal step-size  DT = ', F10.6,', MAG = ', F12.6)
  907 FORMAT(' TAU(0) = ', F12.6,', TAU(1) = ', F12.6, ' in dimensionless units ')
  906 FORMAT(' A0 = ', F12.6, ', A2 = ', F12.6, ' (all in units of Bohr radius)')
  908 FORMAT(' AOSC = ', E12.4)
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
          DO K = 1, NZ
             READ(1,*) (TMP(L), L = 1, 4*SPIN+2)
             SELECT CASE(SPIN)
               CASE(1)
                 PHI(I,J,K,1) = SQRT(TMP(1))*EXP(CI*TMP(2))
                 PHI(I,J,K,2) = SQRT(TMP(3))*EXP(CI*TMP(4))
                 PHI(I,J,K,3) = SQRT(TMP(5))*EXP(CI*TMP(6))
               CASE(2)
                 PHI(I,J,K,1) = SQRT(TMP(1))*EXP(CI*TMP(2))
                 PHI(I,J,K,2) = SQRT(TMP(3))*EXP(CI*TMP(4))
                 PHI(I,J,K,3) = SQRT(TMP(5))*EXP(CI*TMP(6))
                 PHI(I,J,K,4) = SQRT(TMP(7))*EXP(CI*TMP(8))
                 PHI(I,J,K,5) = SQRT(TMP(9))*EXP(CI*TMP(10))
             END SELECT
          END DO
       END DO
    END DO
  END IF
  
  IF(SWITCH_IM.EQ.1) CALL NORMT(PHI, TNORM)

  MS = 2*SPIN+1

  CALL NORMC(PHI, NORM)
  CALL RAD(PHI, RMS)
  CALL ENERGY(PHI, MU, EN, MZ)
   
  SELECT CASE(SPIN)
    CASE(1)
      WRITE (2, 1001)
      WRITE (2, 1002)
      WRITE (2, 1001)
      WRITE (2, 1003) SUM(NORM), EN/2.0D0, (MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,NY2,NZ2,K)), K = 1, MS)
    CASE(2)
      WRITE (2, 1011)
      WRITE (2, 1012)
      WRITE (2, 1011)
      WRITE (2, 1013) SUM(NORM), EN/2.0D0, (MU(K)/2.0D0, K = 1, MS), (ABS(PHI(NX2,NY2,NZ2,K)), K = 1, MS)
  END SELECT

  1001 FORMAT (23X, 74('-'))
  1002 FORMAT (23X, 'Norm', 7X,'EN', 7X, 'MU1', 7X, 'MU2', 7X, 'MU3', 6X,&
                    'phi1', 6X, 'phi2', 6X,'phi3')
  1003 FORMAT ('Initial: ', 8X, 1F11.4, 4F10.5, 3F10.5)
  1011 FORMAT (22X, 117('-'))
  1012 FORMAT (23X, 'Norm', 7X,'EN', 7X, 'MU1', 7X, 'MU2', 7X, 'MU3', 7X, 'MU4', 7X, 'MU5', 7X,&
                    'phi1', 7X, 'phi2', 5X,'phi3',5X, 'phi4', 6X,'phi5')
  1013 FORMAT ('Initial : ', 7X, 1F11.4, 6F10.5, 5F10.5)

  !Initializing SIG(1) and SIG(2)
  SIG(1) = 20.0D0
  SIG(2) = 20.0D0

  CALL CREATE_PLANS()

  WRITE(3,'(F14.6,6F18.10)') 2.0D0*II*DT, EN/2.0D0, (RMS(L), L = 1, MS)
  WRITE(4,'(F14.6,7F18.10)') 2.0D0*II*DT, (NORM(L), L = 1, MS), SUM(NORM), MZ

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
             SIGMA(2) = SQRT(1.0D0 - MAG**2)/SQRT(NORM(2)+SQRT(4.0D0*&
                      (1.0D0-MAG**2)*NORM(1)*NORM(3)+MAG**2*NORM(2)**2)+EPS)
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
    
         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L)
         DO L = 1, MS
            DO K = 1, NZ
               DO J = 1, NY
                  DO I = 1, NX
                     PHI(I,J,K,L) = SIGMA(L)*PHI(I,J,K,L)
                  END DO
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
       ! In tmp_solution_file.dat, density of first component followed by its                      !
       ! phase is written in first and second columns, respectively, then                          !
       ! same for the rest of the components.                                                      !
       !                                                                                           !
       ! In convergence.dat, time, wavefunction convergence and energy convergence                 ! 
       ! parameter are written in  first and second columns, respectively.                         !
       !-------------------------------------------------------------------------------------------!

        WRITE(3,'(F14.6,6F18.10)') 2.0D0*II*DT, EN/2.0D0, (RMS(L), L = 1, MS)
        WRITE(4,'(F14.6,7F18.10)') 2.0D0*II*DT, (NORM(L), L = 1, MS), SUM(NORM), MZ
     END IF
     EN_OLD = EN
     PHI_OLD(:,:,:,:) = PHI(:,:,:,:)

     IF(MOD(II,NSTP).EQ.0)THEN
       WRITE (2, 1005) SUM(NORM), EN/2.0D0, (MU(L)/2.0D0, L = 1, MS), (ABS(PHI(NX2,NY2,NZ2,L)), L = 1, MS)
       OPEN (7, FILE='./OUTPUT/tmp_solution_file_spin'//TRIM(ADJUSTL(SPIN_ID))//'.dat',& 
                STATUS='replace',FORM='FORMATTED',ACTION='WRITE')
       DO I = 1, NX
          DO J = 1, NY
             DO K = 1, NZ
                WRITE(7,'(10F16.10)')  (ABS(PHI(I,J,K,L))**2, ATAN2(AIMAG(PHI(I,J,K,L)),REAL(PHI(I,J,K,L))), L = 1, MS)
             END DO
          END DO
       END DO
       CLOSE(7) 

       OPEN (11, FILE='./OUTPUT/tmp_solution_file_spin'//TRIM(ADJUSTL(SPIN_ID))//'_xy.dat',& 
                 STATUS='replace',FORM='FORMATTED',ACTION='WRITE')
       DO I = 1, NX
          DO J = 1, NY
             WRITE(11,'(2F10.4,10F14.8)') X(I), Y(J), (ABS(PHI(I,J,NZ2,L))**2,&
                                  ATAN2(AIMAG(PHI(I,J,NZ2,L)),REAL(PHI(I,J,NZ2,L))),L=1,MS)     
          END DO
          WRITE(11,'(2F10.4,10F14.8)')
       END DO
       CLOSE(11)

       OPEN (12, FILE='./OUTPUT/tmp_solution_file_spin'//TRIM(ADJUSTL(SPIN_ID))//'_xz.dat',& 
                 STATUS='replace',FORM='FORMATTED',ACTION='WRITE')
       DO I = 1, NX
          DO J = 1, NZ
             WRITE(12,'(2F10.4,10F14.8)') X(I), Z(J), (ABS(PHI(I,NY2,J,L))**2,&
                                    ATAN2(AIMAG(PHI(I,NY2,J,L)),REAL(PHI(I,NY2,J,L))),L=1,MS)
          END DO
          WRITE(12,'(2F10.4,10F14.8)')
       END DO
       CLOSE(12)
     END IF
  END DO !! Iteration loop

  CALL DESTROY_PLANS()
  CALL ENERGY(PHI, MU, EN, MZ)
  CALL RAD(PHI, RMS)

  IF (SPIN.EQ.1) THEN
      WRITE (2, 1001)
  ELSE
     WRITE (2, 1011)
  END IF

  WRITE (2, 1006) II-1, SUM(NORM), EN/2.0D0,(MU(L)/2.0D0, L = 1, MS), (ABS(PHI(NX2,NY2,NZ2,L)), L = 1, MS)
  1005 FORMAT('After NSTP iter.:    ', 1F7.4, 6F10.5, 5F10.5)
  1006 FORMAT('After ', I8 ,' iter.:', 1F7.4, 6F10.5, 5F10.5)

  IF (SPIN.EQ.1) THEN
     WRITE (2, 1001)
  ELSE
     WRITE (2, 1011)
  END IF
  WRITE(2,*)

  !-------------------------------------------------------------------------------------------------
  ! In solution_file.dat, density of first component followed by its                               ! 
  ! phase is written in first and second columns, respectively, followed                           !
  ! by the same for rest of the components.                                                        !
  !-------------------------------------------------------------------------------------------------

  DO I = 1, NX
     DO J = 1, NY
        DO K = 1, NZ
           WRITE(8,'(10F16.10)') (ABS(PHI(I,J,K,L))**2, ATAN2(AIMAG(PHI(I,J,K,L)),REAL(PHI(I,J,K,L))), L = 1, MS)
        END DO
     END DO
  END DO
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
END PROGRAM CGPE3D

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   Dynamic Memory Allocation to the variables and arrays                          !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE ALLOCATE_MEM()
  USE DOUBLE_PRECISION
  USE BASIC_DATA, ONLY : NX, NY, NZ
  USE CGPE_DATA, ONLY : SPIN, X, X2, Y, Y2, Z, Z2, V, R2, KX,&
                       KY, KZ, PHI, PHIF, PHI_OLD, SX, SY, SZ

  USE FFTW_DATA, ONLY : FFTFXYZ, FFTBXYZ

  IMPLICIT NONE
  INTEGER :: MS, STAT

  MS = 2*SPIN+1
  ALLOCATE(X(1:NX),X2(1:NX),Y(1:NY), Y2(1:NY), Z(1:NZ), Z2(1:NZ), V(1:NX,1:NY,1:NZ), R2(1:NX,1:NY,1:NZ), &
         KX(1:NX), KY(1:NY),  KZ(1:NZ), STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(PHI(1:NX,1:NY,1:NZ,1:MS),PHIF(1:NX,1:NY,1:NZ,1:MS), PHI_OLD(1:NX,1:NY,1:NZ,1:MS), STAT = STAT)
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(SX(1:MS,1:MS),SY(1:MS,1:MS),SZ(1:MS))
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
  ALLOCATE(FFTFXYZ(NX,NY,NZ), FFTBXYZ(NX,NY,NZ))
  IF(STAT.NE.0) WRITE(100,*) 'Error in memory allocation.'
END SUBROUTINE ALLOCATE_MEM

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                       Memory Deallocation of the variables and arrays                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DEALLOCATE_MEM()
  USE DOUBLE_PRECISION
  USE CGPE_DATA, ONLY : X, X2, Y,Y2,Z, Z2, V, R2, KX, KY, KZ, PHI, PHIF, PHI_OLD, SX, SY, SZ
  USE FFTW_DATA, ONLY : FFTFXYZ, FFTBXYZ

  IMPLICIT NONE
  INTEGER :: STAT

  DEALLOCATE(X, X2, Y, Y2, Z, Z2, V, R2, KX, KY, KZ, STAT = STAT)
  DEALLOCATE(PHI, PHIF, PHI_OLD, STAT = STAT)
  DEALLOCATE(SX,SY,SZ)
  DEALLOCATE(FFTFXYZ, FFTBXYZ)
END SUBROUTINE DEALLOCATE_MEM

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Defines the SX, SY, SZ, spatial and Fourier grids, trapping potential,              !
!              and initializes the component wave functions                                        !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INITIALIZE()
  USE BASIC_DATA
  USE CGPE_DATA
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ):: TMP3D, TFDENSITY
  INTEGER :: I, J, K
  
  SELECT CASE (SPIN)
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
  DO I = 1, 1+NZ/2
     KZ(I) = (I-1.0d0)*2.0D0*PI/LZ
  END DO

  DO I = 1,NZ/2 -1
     KZ(I+1+NZ/2)=-KZ(1-I+NZ/2)
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
  
  !$OMP SECTION
  DO I = 1, NZ
     Z(I) = (-1.0d0 +2.0D0*(I-1.0D0)/NZ)*LZ/2.0D0
     Z2(I) = Z(I)*Z(I)
  END DO
  !$OMP END PARALLEL SECTIONS 

  !$OMP PARALLEL DO PRIVATE(I,J,K) 
  DO K = 1, NZ
     DO J = 1, NY
        DO I = 1, NX
           R2(I,J,K) = X2(I) + Y2(J)+ Z2(K)
           V(I,J,K) = ALPHAX*ALPHAX*X2(I) + ALPHAY*ALPHAY*Y2(J) + ALPHAZ*ALPHAZ*Z2(K)
           TMP3D(I,J,K) = ALPHAX*X2(I)+ ALPHAY*Y2(J) + ALPHAZ*Z2(K)
           SELECT CASE(SPIN)
 
             CASE(1)
               SELECT CASE (OPTION_FPC)
                 CASE(1) !Initial guess wavefunctions for Ferromagnetic interactions
                   IF(((15.0D0*(TAU(0)+TAU(1))/(8.0D0*PI))**(0.4D0)- ALPHAX*ALPHAX*X2(I)&
                     -ALPHAY*ALPHAY*Y2(J)- ALPHAZ*ALPHAZ*Z2(K)).GE.0.0D0)THEN
                      TFDENSITY(I,J,K) = ((15.0D0*(TAU(0)+TAU(1))/(8.0D0*PI))**(0.4D0)- ALPHAX*ALPHAX*X2(I) -&
                                          ALPHAY*ALPHAY*Y2(J) - ALPHAZ*ALPHAZ*Z2(K))/((TAU(0)+TAU(1)))
                   ELSE
                     TFDENSITY(I,J,K) = 0.0D0
                   END IF

                   IF(((15.0D0*(TAU(0)+TAU(1))/(8.0D0*PI))**(0.2D0)).GE.10.0D0)THEN
                     PHI(I,J,K,1) = ((1.0D0+MAG)/2.0D0)*SQRT(TFDENSITY(I,J,K))
                     PHI(I,J,K,2) = SQRT((1.0D0-MAG**2)/2.0D0)*SQRT(TFDENSITY(I,J,K))
                     PHI(I,J,K,3) = ((1.0D0-MAG)/2.0D0)*SQRT(TFDENSITY(I,J,K))
                   ELSE
                     PHI(I,J,K,1) = ((1.0D0+MAG)/2.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                     PHI(I,J,K,2) = SQRT((1.0D0-MAG**2)/2.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                     PHI(I,J,K,3) = ((1.0D0-MAG)/2.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   END IF

                 CASE(2) !Intitial guess wavefunctions for antiferromagnetic interactions 
                   PHI(I,J,K,1) = (SQRT((1.0D0+MAG)/2.0D0))*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,2) = 0.0D0
                   PHI(I,J,K,3) = (SQRT((1.0D0-MAG)/2.0D0))*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))

                 CASE DEFAULT!For choosing Gaussian intial guess wavefunction
                   PHI(I,J,K,1) = (EXP(-(X2(I)+Y2(J)+Z2(K))/2.0D0)/(PI**0.75D0))
                   PHI(I,J,K,2) = (EXP(-(X2(I)+Y2(J)+Z2(K))/2.0D0)/(PI**0.75D0))
                   PHI(I,J,K,3) = (EXP(-(X2(I)+Y2(J)+Z2(K))/2.0D0)/(PI**0.75D0))
               END SELECT
             CASE(2)
               SELECT CASE (OPTION_FPC)
                 CASE(1)!Initial guess wavefunctions for Ferromagnetic interactions
                   IF(((15.0D0*(TAU(0)+4.0D0*TAU(1))/(8.0D0*PI))**(0.4D0)- ALPHAX*ALPHAX*X2(I) &
                     - ALPHAY*ALPHAY*Y2(J)- ALPHAZ*ALPHAZ*Z2(K)).GE.0.0D0)THEN
                       TFDENSITY(I,J,K) = ((15.0D0*(TAU(0)+4.0D0*TAU(1))/(8.0D0*PI))**(0.4D0)- ALPHAX*ALPHAX*X2(I) -&
                                           ALPHAY*ALPHAY*Y2(J) - ALPHAZ*ALPHAZ*Z2(K))/((TAU(0)+4.0D0*TAU(1)))
                   ELSE
                     TFDENSITY(I,J,K) = 0.0D0
                   END IF
       
                   IF(((4.0D0*(TAU(0)+4.0D0*TAU(1))/(2.0D0*PI))**(0.25D0)).GE.10.0D0)THEN
                     PHI(I,J,K,1) = (((2.0D0+MAG)**2)/16.0D0)*SQRT(TFDENSITY(I,J,K))
                     PHI(I,J,K,2) = ((SQRT(4.0D0-MAG**2)*(2.0D0+MAG))/8.0D0)*SQRT(TFDENSITY(I,J,K))
                     PHI(I,J,K,3) = ((SQRT(1.5D0)*(4.0D0-MAG**2))/8.0D0)*SQRT(TFDENSITY(I,J,K))
                     PHI(I,J,K,4) = ((SQRT(4.0D0-MAG**2)*(2.0D0-MAG))/8.0D0)*SQRT(TFDENSITY(I,J,K))
                     PHI(I,J,K,5) = (((2.0D0-MAG)**2)/16.0D0)*SQRT(TFDENSITY(I,J,K))
                   ELSE
                     PHI(I,J,K,1) = (((2.0D0+MAG)**2)/16.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                     PHI(I,J,K,2) = ((SQRT(4.0D0-MAG**2)*(2.0D0+MAG))/8.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                     PHI(I,J,K,3) = ((SQRT(1.5D0)*(4.0D0-MAG**2))/8.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                     PHI(I,J,K,4) = ((SQRT(4.0D0-MAG**2)*(2.0D0-MAG))/8.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                     PHI(I,J,K,5) = (((2.0D0-MAG)**2)/16.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   END IF

                 CASE(2)!Intitial guess wavefunctions for polar interactions 
                   PHI(I,J,K,1) = (SQRT(2.0D0+MAG)/2.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,2) = 0.0D0
                   PHI(I,J,K,3) = 0.0D0
                   PHI(I,J,K,4) = 0.0D0
                   PHI(I,J,K,5) = (SQRT(2.0D0-MAG)/2.0D0)*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))

                 CASE(3)!Intitial guess wavefunctions for cyclic interactions 
                   PHI(I,J,K,1) = (SQRT((1.0D0+MAG)/3.0D0))*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,2) = 0.0D0
                   PHI(I,J,K,3) = 0.0D0
                   PHI(I,J,K,4) = (SQRT((2.0D0-MAG)/3.0D0))*EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,5) = 0.0D0

                 CASE DEFAULT !For choosing Gaussian intial guess wavefunction
                   PHI(I,J,K,1) = EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,2) = EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,3) = EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,4) = EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
                   PHI(I,J,K,5) = EXP(-TMP3D(I,J,K)/2.0D0)/((PI)**(0.75D0))
               END SELECT
           END SELECT
        END DO
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
  PLANFXYZ = FFTW_PLAN_DFT_3D(NZ,NY,NX,FFTFXYZ,FFTBXYZ,FFTW_FORWARD,FFTW_ESTIMATE)
  PLANBXYZ = FFTW_PLAN_DFT_3D(NZ,NY,NX,FFTBXYZ,FFTFXYZ,FFTW_BACKWARD,FFTW_ESTIMATE)
END SUBROUTINE CREATE_PLANS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                    Destroys FFTW plans for forward and backward transforms                       !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE DESTROY_PLANS()
  USE FFTW3
  USE FFTW_DATA
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  CALL FFTW_DESTROY_PLAN(PLANFXYZ)
  CALL FFTW_DESTROY_PLAN(PLANBXYZ)
  CALL FFTW_CLEANUP_THREADS()
END SUBROUTINE DESTROY_PLANS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    Solves the Hamiltonian corresponding to Kinetic Energy and diagonal spin orbit coupling terms !
!    in Fourier space                                                                              !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KE()
  USE BASIC_DATA, ONLY : CI, NX, NY, NZ,CDT
  USE CGPE_DATA, ONLY :  PHIF, KX, KY, KZ, SPIN
  USE SOC_DATA, ONLY : GAMMAZ
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, L, MS
  MS = 2*SPIN+1

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L)
  DO L = -SPIN, SPIN, 1
     DO K = 1, NZ
        DO J = 1, NY
           DO I = 1, NX
              PHIF(I,J,K,L+SPIN+1) = EXP(-CI*CDT*(KX(I)*KX(I) + KY(J)*KY(J)+KZ(K)*KZ(K)&
                                    -2.0D0*L*GAMMAZ*KZ(K)))*PHIF(I,J,K,L+SPIN+1) 
           END DO
        END DO
     END DO
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE KE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       Solves the Hamiltonian corresponding to spin-orbit coupling terms in Fourier space         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SOC()
  USE BASIC_DATA, ONLY : CI, NX, NY, NZ, EPS, CDT
  USE CGPE_DATA, ONLY : KX, KY, PHIF, SPIN
  USE SOC_DATA, ONLY : GAMMAX, GAMMAY
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(1:2*SPIN+1,1:2*SPIN+1) :: A
  REAL(KIND=DBL) :: ALPHA, BETA
  INTEGER :: I, J, K, L, MS

  MS = 2*SPIN+1
  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L,ALPHA,BETA,A)
      DO I = 1, NX
         DO J = 1, NY
            ALPHA= 2.0D0*SQRT((GAMMAX**2)*KX(I)**2+(GAMMAY**2)*KY(J)**2)
            IF(ALPHA.LE.EPS) THEN
              BETA = 0.0D0
            ELSE
              BETA= ATAN2(GAMMAY*KY(J),GAMMAX*KX(I))
            END IF
    
            A(1,1) = 0.5D0*EXP(-2.0D0*CI*BETA)
            A(2,1) = -SQRT(0.5D0)*EXP(-CI*BETA)
            A(3,1) = 0.5D0

            A(1,2) = -SQRT(0.5D0)*EXP(-2.0D0*CI*BETA)
            A(2,2) = 0.0D0
            A(3,2) = SQRT(0.5D0)

            A(1,3) = A(1,1)
            A(2,3) = -A(2,1)
            A(3,3) = A(3,1)
          
            DO K = 1, NZ
               PHIF(I,J,K,1:MS) = MATMUL(CONJG(TRANSPOSE(A(1:MS,1:MS))), PHIF(I,J,K,1:MS))
               DO L = -SPIN, SPIN, 1
                  PHIF(I,J,K,L+SPIN+1) = EXP(-CI*CDT*(L*ALPHA)) * PHIF(I,J,K,L+SPIN+1)
               END DO
               PHIF(I,J,K,1:MS) = MATMUL(A(1:MS,1:MS), PHIF(I,J,K,1:MS))
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L,ALPHA,BETA,A)
      DO I = 1, NX
         DO J = 1, NY   
            ALPHA= 2.0D0*SQRT((GAMMAX**2)*KX(I)**2+(GAMMAY**2)*KY(J)**2)
            IF(ALPHA.LE.EPS) THEN
              BETA = 0.0D0
            ELSE 
              BETA= ATAN2(GAMMAY*KY(J),GAMMAX*KX(I))
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
            A(3,3) = -CMPLX(0.5D0,0.0D0,KIND=DBL)!
            A(4,3) = A(3,2)
            A(5,3) = CONJG(A(1,3))

            A(1,4) = -0.5D0*EXP(-2.0D0*CI*BETA)
            A(2,4) = -0.5D0*EXP(-CI*BETA)
            A(3,4) = A(3,2)
            A(4,4) = -CONJG(A(2,4))
            A(5,4) = -CONJG(A(1,4))

            A(1,5) = 0.25D0*EXP(-2.0D0*CI*BETA)
            A(2,5) = 0.5D0*EXP(-CI*BETA)
            A(3,5) = SQRT(0.375D0)
            A(4,5) = CONJG(A(2,5))
            A(5,5) = CONJG(A(1,5))            

            DO K = 1, NZ
               PHIF(I,J,K,1:MS) = MATMUL(CONJG(TRANSPOSE(A(1:MS,1:MS))), PHIF(I,J,K,1:MS))      
               DO L = -SPIN, SPIN, 1
                  PHIF(I,J,K,L+SPIN+1) = EXP(-CI*CDT*(L*ALPHA)) * PHIF(I,J,K,L+SPIN+1)
               END DO
               PHIF(I,J,K,1:MS) = MATMUL(A(1:MS,1:MS), PHIF(I,J,K,1:MS))
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
  END SELECT
END SUBROUTINE SOC

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    Solves the Hamiltonian corresponding to trapping potential and diagonal interaction terms     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SP()
  USE BASIC_DATA, ONLY : CI, NX, NY, NZ, CDT
  USE CGPE_DATA, ONLY :  V, TAU, PHI, SPIN !, C0, C2
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  REAL(KIND=DBL), DIMENSION(1:2*SPIN+1) :: TMP
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ) :: RHO, FZ
  INTEGER :: I, J, K, L, MS

  MS = 2*SPIN+1
  SELECT CASE (SPIN)
    CASE(1)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,TMP,L)
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               RHO(I,J,K) = SUM(ABS(PHI(I,J,K,1:MS))**2)
               TMP(1) = V(I,J,K) + TAU(0)*RHO(I,J,K) + TAU(1)*(ABS(PHI(I,J,K,1))**2 + ABS(PHI(I,J,K,2))**2&
                      -ABS(PHI(I,J,K,3))**2)
               TMP(2) = V(I,J,K) + TAU(0)*RHO(I,J,K) + TAU(1)*(ABS(PHI(I,J,K,1))**2 + ABS(PHI(I,J,K,3))**2)
               TMP(3) = V(I,J,K) + TAU(0)*RHO(I,J,K) + TAU(1)*(ABS(PHI(I,J,K,3))**2 + ABS(PHI(I,J,K,2))**2&
                        -ABS(PHI(I,J,K,1))**2)
               DO L = 1, MS
                  PHI(I,J,K,L) = PHI(I,J,K,L)*EXP((-CI*CDT)*TMP(L))
               END DO
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,TMP,L)
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               RHO(I,J,K) = SUM(ABS(PHI(I,J,K,1:MS))**2)
               FZ(I,J,K)  = 2.0D0*ABS(PHI(I,J,K,1))**2 + ABS(PHI(I,J,K,2))**2&
                            - ABS(PHI(I,J,K,4))**2 - 2.0D0*ABS(PHI(I,J,K,5))**2
               TMP(1) = (V(I,J,K) + TAU(0)*RHO(I,J,K) + 2.0D0*TAU(1)*FZ(I,J,K)&
                        + 0.4D0*TAU(2)*ABS(PHI(I,J,K,5))**2)
               TMP(2) = (V(I,J,K) + TAU(0)*RHO(I,J,K) + TAU(1)*FZ(I,J,K)&
                        + 0.4D0*TAU(2)*ABS(PHI(I,J,K,4))**2) 
               TMP(3) = (V(I,J,K) + TAU(0)*RHO(I,J,K)+ 0.2D0*TAU(2)*ABS(PHI(I,J,K,3))**2)
               TMP(4) = (V(I,J,K) + TAU(0)*RHO(I,J,K) - TAU(1)*FZ(I,J,K) + 0.4D0*TAU(2)*ABS(PHI(I,J,K,2))**2)
               TMP(5) = (V(I,J,K)+TAU(0)*RHO(I,J,K) - 2.0D0*TAU(1)*FZ(I,J,K)&
                        + 0.4D0*TAU(2)*ABS(PHI(I,J,K,1))**2)
               DO L = 1, MS
                  PHI(I,J,K,L) = PHI(I,J,K,L)*EXP((-CI*CDT)*TMP(L))
               END DO           
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
  USE BASIC_DATA, ONLY : NX, NY, NZ
  USE CGPE_DATA, ONLY : SPIN, SX, SY, SZ 
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE
  INTEGER :: I,J,K,MS
  COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT) :: FX, FY, FZ
  MS = 2*SPIN+1
  
  !$OMP PARALLEL DO PRIVATE(I,J,K)
  DO K = 1, NZ
     DO J = 1, NY
        DO I = 1, NX
            FX(I,J,K) = REAL(DOT_PRODUCT(PHI(I,J,K,1:MS), MATMUL(SX(1:MS,1:MS),PHI(I,J,K,1:MS))))
            FY(I,J,K) = REAL(DOT_PRODUCT(PHI(I,J,K,1:MS), MATMUL(SY(1:MS,1:MS),PHI(I,J,K,1:MS))))
            FZ(I,J,K) = DOT_PRODUCT(SZ(1:MS),ABS(PHI(I,J,K,1:MS))**2)
        END DO
     END DO
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE FXYZ

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            Calculates C12, C13, C23, C34, C35 and C45, i.e. the elements of Hamiltonian          !
!            corresponding to off-diagonal interaction terms in spin-2 systems                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
  USE BASIC_DATA, ONLY : NX, NY,NZ
  USE CGPE_DATA, ONLY : TAU
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ,1:5), INTENT(IN) :: PHI
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(IN) :: FMINUS
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
  
  !$OMP PARALLEL WORKSHARE
  C12(:,:,:) = TAU(1)*FMINUS(:,:,:) - 0.4D0*TAU(2)*PHI(:,:,:,4)*CONJG(PHI(:,:,:,5)) 

  C13(:,:,:) = 0.2D0*TAU(2)*PHI(:,:,:,3)*CONJG(PHI(:,:,:,5))

  C23(:,:,:) = (SQRT(1.5D0))*TAU(1)*FMINUS(:,:,:)-0.2D0*TAU(2)*PHI(:,:,:,3)*CONJG(PHI(:,:,:,4)) 

  C34(:,:,:) = (SQRT(1.5D0))*TAU(1)*FMINUS(:,:,:)-0.2D0*TAU(2)*PHI(:,:,:,2)*CONJG(PHI(:,:,:,3)) 

  C35(:,:,:) = 0.2D0*TAU(2)*PHI(:,:,:,1)*CONJG(PHI(:,:,:,3))

  C45(:,:,:) = TAU(1)*FMINUS(:,:,:) - 0.4D0*TAU(2)*PHI(:,:,:,1)*CONJG(PHI(:,:,:,2)) 
  !$OMP END PARALLEL WORKSHARE
END SUBROUTINE MAT_C

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          Calculates the norm of individual components                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NORMC(PHI, NORM)
  USE BASIC_DATA, ONLY : NX, NY, NZ, DX, DY, DZ
  USE CGPE_DATA, ONLY : SPIN
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
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
  
  INTEGER :: I, J, K, L, MS
  REAL(KIND=DBL), DIMENSION(1:NX) :: P
  REAL(KIND=DBL), DIMENSION(1:NY) :: TMP2D
  REAL(KIND=DBL), DIMENSION(1:NZ,1:2*SPIN+1) :: TMP1D
  MS = 2*SPIN+1  

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L,P,TMP2D)
  DO L = 1, MS
     DO K = 1, NZ
        DO J = 1, NY
           DO I = 1, NX
              P(I) = ABS(PHI(I,J,K,L))**2
           END DO
           TMP2D(J) = SIMPSON(P(:), DX)
        END DO
        TMP1D(K,L) = SIMPSON(TMP2D(:), DY)
     END DO
  END DO
  !$OMP END PARALLEL DO

  DO L = 1, MS
     NORM(L) = SIMPSON(TMP1D(:,L), DZ)
  END DO
END SUBROUTINE NORMC

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Calculates the total norm and normalizes the total density to 1                     !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE NORMT(PHI, TNORM)
  USE BASIC_DATA, ONLY : NX, DX, NY, DY, NZ, DZ
  USE CGPE_DATA, ONLY : SPIN
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(INOUT) :: PHI
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

  INTEGER :: I, J, K, MS
  REAL(KIND=DBL), DIMENSION(1:NX) :: P
  REAL(KIND=DBL), DIMENSION(1:NY) :: TMP2D
  REAL(KIND=DBL), DIMENSION(1:NZ) :: TMP1D
  MS = 2*SPIN+1

  !$OMP PARALLEL DO PRIVATE(I,J,K,P,TMP2D)
  DO K = 1, NZ
     DO J = 1, NY
        DO I = 1, NX
           P(I) =  SUM(ABS(PHI(I,J,K,1:MS))**2)
        END DO
        TMP2D(J) = SIMPSON(P(:), DX)
     END DO
     TMP1D(K) = SIMPSON(TMP2D(:), DY)
  END DO
  !$OMP END PARALLEL DO

  TNORM = SQRT(SIMPSON(TMP1D, DZ))
  PHI = PHI/TNORM
END SUBROUTINE NORMT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!              Calculates the root mean square sizes of the (2*SPIN+1) components                  !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE RAD(PHI, RMS)
  USE BASIC_DATA, ONLY : NX, NY, NZ, DX, DY, DZ
  USE CGPE_DATA, ONLY : R2, SPIN
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
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

  INTEGER :: I, J, K, L, MS
  REAL(KIND=DBL), DIMENSION(1:NX) :: P
  REAL(KIND=DBL), DIMENSION(1:NY) :: TMP2D
  REAL(KIND=DBL), DIMENSION(1:NZ,1:2*SPIN+1) :: TMP1D
  MS = 2*SPIN+1

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L,P,TMP2D)
  DO L = 1,  MS
     DO K = 1, NZ
        DO J = 1, NY
           DO I = 1, NX
              P(I) = R2(I,J,K)*ABS(PHI(I,J,K,L))*ABS(PHI(I,J,K,L))
           END DO
           TMP2D(J) = SIMPSON(P(:), DX)
        END DO
        TMP1D(K,L) = SIMPSON(TMP2D(:), DY)
     END DO
  END DO
  !$OMP END PARALLEL DO

  DO L = 1, MS
     RMS(L) = SQRT(SIMPSON(TMP1D(:,L), DZ))
  END DO
END SUBROUTINE RAD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       Calculates the (2*SPIN+1) component chemical potentials, energy, and magnetization         !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE ENERGY(PHI, MU, EN, MZ)
  USE BASIC_DATA, ONLY : NX, NY, NZ, DX, DY, DZ, CI, EPS
  USE CGPE_DATA, ONLY : V, TAU, SPIN
  USE SOC_DATA
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
  REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: MU
  REAL(KIND=DBL), INTENT(OUT) :: EN, MZ

  INTERFACE
    PURE FUNCTION SIMPSON(F, DX)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: F
      REAL(KIND=DBL), INTENT(IN) :: DX
      REAL(KIND=DBL) :: SIMPSON
    END FUNCTION SIMPSON

    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      REAL(KIND=DBL), DIMENSION(1:), INTENT(IN) :: P
      REAL(KIND=DBL), INTENT(IN) :: DX
      REAL(KIND=DBL), DIMENSION(1:SIZE(P)) :: DP
    END FUNCTION DIFF
 
    SUBROUTINE FXYZ(PHI, FX, FY, FZ)
      USE BASIC_DATA, ONLY : NX, NY, NZ,CI
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT) :: FX, FY, FZ
    END SUBROUTINE FXYZ

    SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      USE BASIC_DATA, ONLY : NX, NY,NZ
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ,1:5), INTENT(IN) :: PHI
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(IN) ::  FMINUS
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
    END SUBROUTINE MAT_C

    SUBROUTINE NORMC(PHI, NORM)
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(:), INTENT(OUT) :: NORM
    END  SUBROUTINE NORMC
  END INTERFACE

  INTEGER :: I, J, K, L
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ,1:2*SPIN+1) :: DPX, DPY, DPZ
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ,1:2*SPIN+1) :: DP2, DPX2, DPY2
  REAL(KIND=DBL), DIMENSION(1:NX, 1:NY,1:NZ) :: FX, FY, FZ, RHO
  COMPLEX(KIND=DBL), DIMENSION(1:NX, 1:NY,1:NZ) :: THETA, FMINUS, C12, C13, C23, C34, C35, C45
  REAL(KIND=DBL), DIMENSION(1:NX, 2*SPIN+3) :: TMP3D
  REAL(KIND=DBL), DIMENSION(1:NY, 2*SPIN+3) :: TMP2D
  REAL(KIND=DBL), DIMENSION(1:2*SPIN+1) :: NORM
  REAL(KIND=DBL), DIMENSION(1:NZ,1:2*SPIN+1+2) :: TMP1D

  CALL FXYZ(PHI, FX, FY, FZ)

  CALL NORMC(PHI, NORM)
 
  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L)
  DO L = 1, 2*SPIN+1
     DO K = 1, NZ
        DO J = 1, NY
           DPX2(:,J,K,L) = DIFF(REAL(PHI(:,J,K,L)), DX)**2 + DIFF(AIMAG(PHI(:,J,K,L)), DX)**2
        END DO
        DO I = 1, NX
           DPY2(I,:,K,L) = DIFF(REAL(PHI(I,:,K,L)), DY)**2 + DIFF(AIMAG(PHI(I,:,K,L)), DY)**2
        END DO
     END DO
  END DO
  !$OMP END PARALLEL DO   
  
  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,L)
  DO L = 1, 2*SPIN+1
     DO J = 1, NY
         DO I = 1, NX
            DP2(I,J,:,L) = DPX2(I,J,:,L)+DPY2(I,J,:,L)+&
                           DIFF(REAL(PHI(I,J,:,L)), DZ)**2 + DIFF(AIMAG(PHI(I,J,:,L)), DZ)**2

            DPZ(I,J,:,L) = DIFF(REAL(PHI(I,J,:,L)), DZ) + CI*DIFF(AIMAG(PHI(I,J,:,L)), DZ)
         END DO
     END DO
  END DO
  !$OMP END PARALLEL DO 

  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(I,J,K,L)
  DO L = 1, 2*SPIN+1
     DO K = 1, NZ
        DO J = 1, NY
           DPX(:,J,K,L) = DIFF(REAL(PHI(:,J,K,L)), DX) + CI*DIFF(AIMAG(PHI(:,J,K,L)), DX)
        END DO
        DO I = 1, NX
           DPY(I,:,K,L) = DIFF(REAL(PHI(I,:,K,L)), DY) + CI*DIFF(AIMAG(PHI(I,:,K,L)), DY)
        END DO
     END DO
  END DO
  !$OMP END PARALLEL DO 
  
  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO PRIVATE(I,J,K,TMP3D,TMP2D)
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               RHO(I,J,K) = SUM(ABS(PHI(I,J,K,1:2*SPIN+1))**2)
               TMP3D(I,1) = REAL(V(I,J,K)*RHO(I,J,K) + TAU(0)*RHO(I,J,K)**2/2.0D0 +&
                            TAU(1)*(ABS(FX(I,J,K))**2 + ABS(FY(I,J,K))**2 + ABS(FZ(I,J,K))**2)/2.0D0+&
                            DP2(I,J,K,1) + DP2(I,J,K,2) + DP2(I,J,K,3) +&
                            2.0D0*(SQRT(0.5D0))*(CONJG(PHI(I,J,K,1))*(-CI*GAMMAX*DPX(I,J,K,2)&
                            -GAMMAY*DPY(I,J,K,2))+CONJG(PHI(I,J,K,2))*(-CI*GAMMAX*DPX(I,J,K,3)-GAMMAY*DPY(I,J,K,3))+&
                            CONJG(PHI(I,J,K,2))*(-CI*GAMMAX*DPX(I,J,K,1)+GAMMAY*DPY(I,J,K,1))&
                            +CONJG(PHI(I,J,K,3))*(-CI*GAMMAX*DPX(I,J,K,2)+GAMMAY*DPY(I,J,K,2)))+&
                            2.0D0*(CONJG(PHI(I,J,K,1))*(-CI*GAMMAZ*DPZ(I,J,K,1))&
                            +CONJG(PHI(I,J,K,3))*(CI*GAMMAZ*DPZ(I,J,K,3))))

               TMP3D(I,2) = ABS(PHI(I,J,K,1))**2-ABS(PHI(I,J,K,3))**2

               TMP3D(I,3) = REAL(V(I,J,K)*ABS(PHI(I,J,K,1))*ABS(PHI(I,J,K,1)) + DP2(I,J,K,1) +&
                            TAU(0)*(ABS(PHI(I,J,K,1))**2+ABS(PHI(I,J,K,2))**2+ABS(PHI(I,J,K,3))**2)*ABS(PHI(I,J,K,1))**2+&
                            TAU(1)*(ABS(PHI(I,J,K,1))**2+ABS(PHI(I,J,K,2))**2-ABS(PHI(I,J,K,3))**2)*ABS(PHI(I,J,K,1))**2+&
                            TAU(1)*CONJG(PHI(I,J,K,3))*PHI(I,J,K,2)**2*CONJG(PHI(I,J,K,1))+&
                            SQRT(2.0D0)*(CONJG(PHI(I,J,K,1))*(-CI*GAMMAX*DPX(I,J,K,2)&
                            -GAMMAY*DPY(I,J,K,2)-CI*SQRT(2.0D0)*GAMMAZ*DPZ(I,J,K,1))))

               TMP3D(I,4) = REAL(V(I,J,K)*ABS(PHI(I,J,K,2))*ABS(PHI(I,J,K,2)) + DP2(I,J,K,2) +&
                            TAU(0)*(ABS(PHI(I,J,K,1))**2+ABS(PHI(I,J,K,2))**2+ABS(PHI(I,J,K,3))**2)*ABS(PHI(I,J,K,2))**2+&
                            TAU(1)*(ABS(PHI(I,J,K,1))**2+ABS(PHI(I,J,K,3))**2)*ABS(PHI(I,J,K,2))**2+&
                            2.0D0*TAU(1)*PHI(I,J,K,1)*PHI(I,J,K,3)*CONJG(PHI(I,J,K,2))**2+&
                            SQRT(2.0D0)*(CONJG(PHI(I,J,K,2))*(-CI*GAMMAX*DPX(I,J,K,3)-GAMMAY*DPY(I,J,K,3))+&
                            CONJG(PHI(I,J,K,2))*(-CI*GAMMAX*DPX(I,J,K,1)+GAMMAY*DPY(I,J,K,1))))

               TMP3D(I,5) = REAL(V(I,J,K)*ABS(PHI(I,J,K,3))*ABS(PHI(I,J,K,3)) + DP2(I,J,K,3) +&
                            TAU(0)*(ABS(PHI(I,J,K,1))**2+ABS(PHI(I,J,K,2))**2+ABS(PHI(I,J,K,3))**2)*ABS(PHI(I,J,K,3))**2+&
                            TAU(1)*(ABS(PHI(I,J,K,2))**2+ABS(PHI(I,J,K,3))**2-ABS(PHI(I,J,K,1))**2)*ABS(PHI(I,J,K,3))**2+&
                            TAU(1)*CONJG(PHI(I,J,K,1))*PHI(I,J,K,2)**2*CONJG(PHI(I,J,K,3))+&
                            SQRT(2.0D0)*(CONJG(PHI(I,J,K,3))*(-CI*GAMMAX*DPX(I,J,K,2)&
                            +GAMMAY*DPY(I,J,K,2)+CI*GAMMAZ*SQRT(2.0D0)*DPZ(I,J,K,3))))

            END DO
            DO I = 1, 2*SPIN+3
               TMP2D(J,I) = SIMPSON(TMP3D(:,I), DX)
            END DO
         END DO
         DO I = 1, 2*SPIN+3
            TMP1D(K,I) = SIMPSON(TMP2D(:,I),DY)
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      FMINUS = FX - CI * FY
      CALL MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      !$OMP PARALLEL DO PRIVATE(I,J,K,TMP3D,TMP2D)
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               RHO(I,J,K) = SUM(ABS(PHI(I,J,K,1:2*SPIN+1))**2)
               THETA(I,J,K) = (2.0D0*PHI(I,J,K,1)*PHI(I,J,K,5)-2.0D0*PHI(I,J,K,2)*PHI(I,J,K,4)+PHI(I,J,K,3)**2)/SQRT(5.0D0)
               TMP3D(I,1) = REAL(V(I,J,K)*RHO(I,J,K) + DP2(I,J,K,1) + &
                            DP2(I,J,K,2) + DP2(I,J,K,3) + DP2(I,J,K,4) + DP2(I,J,K,5) + TAU(0)*RHO(I,J,K)**2/2.0D0 +&
                            TAU(1)*(ABS(FX(I,J,K))**2 + ABS(FY(I,J,K))**2 + ABS(FZ(I,J,K))**2)/2.0D0+&
                            TAU(2)*(THETA(I,J,K)*CONJG(THETA(I,J,K)))/2.0D0&
                            -CI*CONJG(PHI(I,J,K,1))*(2.0D0*GAMMAX*DPX(I,J,K,2)& 
                            -CI*2.0D0*GAMMAY*DPY(I,J,K,2)+ 4.0D0*GAMMAZ*DPZ(I,J,K,1)) -&
                            CI*CONJG(PHI(I,J,K,2))*(GAMMAX*(2.0D0*DPX(I,J,K,1)& 
                            +SQRT(6.0D0)*DPX(I,J,K,3)) + CI*GAMMAY*(2.0D0*DPY(I,J,K,1)-&
                            SQRT(6.0D0)*DPY(I,J,K,3))+ 2.0D0*GAMMAZ*DPZ(I,J,K,2)) &
                            -CI*CONJG(PHI(I,J,K,3))*(GAMMAX*SQRT(6.0D0)*(DPX(I,J,K,2) +&
                            DPX(I,J,K,4)) + CI*GAMMAY*SQRT(6.0D0)*(DPY(I,J,K,2) - DPY(I,J,K,4))) -&
                            CI*CONJG(PHI(I,J,K,4))*(GAMMAX*(SQRT(6.0D0)*DPX(I,J,K,3) +2.0D0*DPX(I,J,K,5)) +&
                            CI*GAMMAY*(SQRT(6.0D0)*DPY(I,J,K,3) - 2.0D0*DPY(I,J,K,5)) - 2.0D0*GAMMAZ*DPZ(I,J,K,4)) -&
                            CI*CONJG(PHI(I,J,K,5))*(2.0D0*GAMMAX*DPX(I,J,K,4)&  
                            +CI*2.0D0*GAMMAY*DPY(I,J,K,4) - 4.0D0*GAMMAZ*DPZ(I,J,K,5)))

               TMP3D(I,2) = 2.0D0*ABS(PHI(I,J,K,1))**2 + ABS(PHI(I,J,K,2))**2 - ABS(PHI(I,J,K,4))**2 -&
                            2.0D0*ABS(PHI(I,J,K,5))**2

               TMP3D(I,3) = REAL((V(I,J,K)+TAU(0)*RHO(I,J,K) + 2.0D0*TAU(1)*FZ(I,J,K)+&
                            0.4D0*TAU(2)*ABS(PHI(I,J,K,5))**2)*ABS(PHI(I,J,K,1))**2+&
                            (C12(I,J,K)*PHI(I,J,K,2)+C13(I,J,K)*PHI(I,J,K,3))*CONJG(PHI(I,J,K,1))-&
                            CI*2.0D0*(CONJG(PHI(I,J,K,1))*GAMMAX*DPX(I,J,K,2)-CI*CONJG(PHI(I,J,K,1))*GAMMAY*DPY(I,J,K,2)+&
                            2.0d0*CONJG(PHI(I,J,K,1))*GAMMAZ*DPZ(I,J,K,1)) + DP2(I,J,K,1))
 
               TMP3D(I,4) = REAL((V(I,J,K)+TAU(0)*RHO(I,J,K) + TAU(1)*FZ(I,J,K)+&
                            0.4D0*TAU(2)*ABS(PHI(I,J,K,4))**2)*ABS(PHI(I,J,K,2))**2+&
                            (CONJG(C12(I,J,K))*PHI(I,J,K,1)+C23(I,J,K)*PHI(I,J,K,3))*CONJG(PHI(I,J,K,2))&
                            -CI*2.0D0*(CONJG(PHI(I,J,K,2))*GAMMAX*DPX(I,J,K,1)+CI*CONJG(PHI(I,J,K,2))*GAMMAY*DPY(I,J,K,1)+&
                            SQRT(1.5D0)*(CONJG(PHI(I,J,K,2))*GAMMAX*DPX(I,J,K,3) -CI*CONJG(PHI(I,J,K,2))*GAMMAY*DPY(I,J,K,3))+&
                            CONJG(PHI(I,J,K,2))*GAMMAZ*DPZ(I,J,K,2))+DP2(I,J,K,2))

               TMP3D(I,5) = REAL((V(I,J,K)+TAU(0)*RHO(I,J,K) +&
                            0.2D0*TAU(2)*ABS(PHI(I,J,K,3))**2)*ABS(PHI(I,J,K,3))**2+&
                            (CONJG(C13(I,J,K))*PHI(I,J,K,1)+CONJG(C23(I,J,K))*PHI(I,J,K,2) &
                            +C34(I,J,K)*PHI(I,J,K,4)+C35(I,J,K)*PHI(I,J,K,5))*CONJG(PHI(I,J,K,3))&
                            -CI*2.0D0*(SQRT(1.5D0)*(CONJG(PHI(I,J,K,3))*GAMMAX*DPX(I,J,K,2)& 
                            +CI*CONJG(PHI(I,J,K,3))*GAMMAY*DPY(I,J,K,2)+&
                            CONJG(PHI(I,J,K,3))*GAMMAX*DPX(I,J,K,4)-CI*CONJG(PHI(I,J,K,3))*GAMMAY*DPY(I,J,K,4)))+&
                            DP2(I,J,K,3))

               TMP3D(I,6) = REAL((V(I,J,K)+TAU(0)*RHO(I,J,K) - TAU(1)*FZ(I,J,K)+&
                            0.4D0*TAU(2)*ABS(PHI(I,J,K,2))**2)*ABS(PHI(I,J,K,4))**2+&
                            (CONJG(C34(I,J,K))*PHI(I,J,K,3)+C45(I,J,K)*PHI(I,J,K,5))*CONJG(PHI(I,J,K,4))&
                            -CI*2.0D0*(CONJG(PHI(I,J,K,4))*GAMMAX*DPX(I,J,K,5)- CI*CONJG(PHI(I,J,K,4))*GAMMAY*DPY(I,J,K,5)+&
                            SQRT(1.5D0)*(CONJG(PHI(I,J,K,4))*GAMMAX*DPX(I,J,K,3) + CI*CONJG(PHI(I,J,K,4))*GAMMAY*DPY(I,J,K,3))&
                            -CONJG(PHI(I,J,K,4))*GAMMAZ*DPZ(I,J,K,4))+ DP2(I,J,K,4))

               TMP3D(I,7) = REAL((V(I,J,K)+TAU(0)*RHO(I,J,K) - 2.0D0*TAU(1)*FZ(I,J,K)+& 
                            0.4D0*TAU(2)*ABS(PHI(I,J,K,1))**2)*ABS(PHI(I,J,K,5))**2+&
                            (CONJG(C35(I,J,K))*PHI(I,J,K,3)+CONJG(C45(I,J,K))*PHI(I,J,K,4))*CONJG(PHI(I,J,K,5))&
                            -CI*2.0D0*(CONJG(PHI(I,J,K,5))*GAMMAX*DPX(I,J,K,4)+CI*CONJG(PHI(I,J,K,5))*GAMMAY*DPY(I,J,K,4)-&
                            2.0D0*CONJG(PHI(I,J,K,5))*GAMMAZ*DPZ(I,J,K,5)) +DP2(I,J,K,5))

            END DO
            DO I = 1, 2*SPIN+3
               TMP2D(J,I) = SIMPSON(TMP3D(:,I), DX)
            END DO
         END DO
         DO I = 1, 2*SPIN+3
            TMP1D(K,I) = SIMPSON(TMP2D(:,I),DY)
         END DO
    END DO
    !$OMP END PARALLEL DO
  END SELECT
  
  EN = SIMPSON(TMP1D(:,1),DZ)/SUM(NORM)
  MZ = SIMPSON(TMP1D(:,2),DZ)
  DO L = 1, 2*SPIN+1
     MU(L) = SIMPSON(TMP1D(:,L+2),DZ)/(NORM(L)+ EPS)
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
  !improve the root. Stop if the root converges.

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
     WRITE(100,*) 'Newton-Raphason subroutine could not calculate the roots within', TOL_NR, ' tolerance'
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
!               Solves the Hamiltonian corresponding to off-diagonal interaction terms             !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE SE()                                                                                   
  USE BASIC_DATA, ONLY : CI, NX, NY, NZ, CDT, EPS
  USE CGPE_DATA, ONLY : PHI, SPIN, TAU
  USE OMP_LIB
  USE DOUBLE_PRECISION

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE FXYZ(PHI, FX, FY, FZ)
      USE BASIC_DATA, ONLY : NX, NY, NZ,CI
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(:,:,:,:), INTENT(IN) :: PHI
      REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT) :: FX, FY, FZ
    END SUBROUTINE FXYZ

    SUBROUTINE MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)
      USE BASIC_DATA, ONLY : NX, NY,NZ
      USE DOUBLE_PRECISION
      IMPLICIT NONE
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ,1:5), INTENT(IN) :: PHI
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(IN) :: FMINUS
      COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT) :: C12, C13, C23, C34, C35, C45
    END SUBROUTINE MAT_C
  END INTERFACE
  
  INTEGER  :: I, J, K, L
  COMPLEX(KIND=DBL), DIMENSION(1:NX, 1:NY, 1:NZ) :: A1, B1, OMEGA, DELTA, KAPPA
  COMPLEX(KIND=DBL), DIMENSION(1:3) :: TMP
  REAL(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ) :: FX, FY, FZ, FX1, FY1, FZ1
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ) :: FMINUS,&
                                 C12, C13, C23, C34, C35, C45,&
                                 C121, C131, C231, C341, C351, C451, FMINUS1
  COMPLEX(KIND=DBL), DIMENSION(1:NX,1:NY,1:NZ,1:5) :: PHI1
  
  EXTERNAL :: ZHEEV
  INTEGER, PARAMETER :: LDA = 5
  INTEGER, PARAMETER :: LWMAX = 1000
  INTEGER :: INFO , LWORK
  REAL(KIND=DBL) :: RWORK(13), W(5)
  COMPLEX(KIND=DBL) ::  A(1:5,1:5), WORK(LWMAX)
 
  SELECT CASE(SPIN)
    CASE(1)
      !$OMP PARALLEL DO PRIVATE(I,J,K,TMP)
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               A1(I,J,K) =  (TAU(1) * PHI(I,J,K,2)*  CONJG(PHI(I,J,K,3)))
               B1(I,J,K) =  (TAU(1) * PHI(I,J,K,2)*  CONJG(PHI(I,J,K,1)))
               OMEGA(I,J,K) = CDT*SQRT((ABS(A1(I,J,K)))**2+(ABS(B1(I,J,K)))**2)

               IF(ABS(OMEGA(I,J,K)).GT.EPS)THEN
                 DELTA(I,J,K) = (COS(OMEGA(I,J,K)) -1)/(OMEGA(I,J,K))**2
                 KAPPA(I,J,K) = -CI * SIN(OMEGA(I,J,K))/OMEGA(I,J,K)
               ELSE
                 DELTA(I,J,K) = -0.5D0
                 KAPPA(I,J,K) = -CI
               END IF

               TMP(1) = (DELTA(I,J,K)*CDT**2)*(A1(I,J,K)*CONJG(A1(I,J,K))*PHI(I,J,K,1)+&
                        A1(I,J,K)*CONJG(B1(I,J,K))*PHI(I,J,K,3)) + KAPPA(I,J,K)*CDT*A1(I,J,K)&
                        *PHI(I,J,K,2)

               TMP(2) = (DELTA(I,J,K)*CDT**2)*(A1(I,J,K)*CONJG(A1(I,J,K))+B1(I,J,K)*CONJG(B1(I,J,K)))&
                        *PHI(I,J,K,2) + KAPPA(I,J,K)*CDT*&
                        (CONJG(A1(I,J,K))*PHI(I,J,K,1)+CONJG(B1(I,J,K))*PHI(I,J,K,3))

               TMP(3) = (DELTA(I,J,K)*CDT**2)*(CONJG(A1(I,J,K))*B1(I,J,K)*&
                        PHI(I,J,K,1)+B1(I,J,K)*CONJG(B1(I,J,K))*PHI(I,J,K,3)) + KAPPA(I,J,K)*CDT*&
                        B1(I,J,K)*PHI(I,J,K,2)

               PHI(I,J,K,1) = PHI(I,J,K,1) + TMP(1)
               PHI(I,J,K,2) = PHI(I,J,K,2) + TMP(2)
               PHI(I,J,K,3) = PHI(I,J,K,3) + TMP(3)
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
    CASE(2)
      CALL FXYZ(PHI, FX, FY, FZ)

      FMINUS(:,:,:) = FX(:,:,:) - CI * FY(:,:,:)

      CALL MAT_C(PHI, FMINUS, C12, C13, C23, C34, C35, C45)

      !$OMP PARALLEL WORKSHARE
      PHI1(:,:,:,1) = PHI(:,:,:,1) - (CI*CDT)*(C12(:,:,:)*PHI(:,:,:,2)+C13(:,:,:)*PHI(:,:,:,3))
      PHI1(:,:,:,2) = PHI(:,:,:,2) - (CI*CDT)*(CONJG(C12(:,:,:))*PHI(:,:,:,1)+C23(:,:,:)*PHI(:,:,:,3))
      PHI1(:,:,:,3) = PHI(:,:,:,3) - (CI*CDT)*(CONJG(C13(:,:,:))*PHI(:,:,:,1)&
                      +CONJG(C23(:,:,:))*PHI(:,:,:,2)+C34(:,:,:)*PHI(:,:,:,4)+C35(:,:,:)*PHI(:,:,:,5))
      PHI1(:,:,:,4) = PHI(:,:,:,4) - (CI*CDT)*(CONJG(C34(:,:,:))*PHI(:,:,:,3)+C45(:,:,:)*PHI(:,:,:,5))
      PHI1(:,:,:,5) = PHI(:,:,:,5) - (CI*CDT)*(CONJG(C35(:,:,:))*PHI(:,:,:,3)+CONJG(C45(:,:,:))*PHI(:,:,:,4))
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
      DO K = 1, NZ
         DO J = 1, NY
            DO I = 1, NX
               A(1,2) = C12(I,J,K)
               A(1,3) = C13(I,J,K)
    
               A(2,1) = CONJG(C12(I,J,K))
               A(2,3) = C23(I,J,K)

               A(3,1) = CONJG(C13(I,J,K))
               A(3,2) = CONJG(C23(I,J,K))
               A(3,4) = C34(I,J,K)
               A(3,5) = C35(I,J,K)
 
               A(4,3) = CONJG(C34(I,J,K))
               A(4,5) = C45(I,J,K)

               A(5,3) = CONJG(C35(I,J,K))
               A(5,4) = CONJG(C45(I,J,K))
  
               LWORK = -1
               CALL ZHEEV( 'Vectors', 'Lower', 5, A, LDA, W, WORK, LWORK, RWORK, INFO )
               LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
               CALL ZHEEV( 'Vectors', 'Lower', 5, A, LDA, W, WORK, LWORK, RWORK, INFO )
       
               PHI(I,J,K,1:5) = MATMUL(CONJG(TRANSPOSE(A(1:5,1:5))), PHI(I,J,K,1:5))
               DO L = 1, 5
                  PHI(I,J,K,L) = EXP(-CI*CDT*W(L)) * PHI(I,J,K,L)
               END DO
               PHI(I,J,K,1:5) = MATMUL(A(1:5,1:5), PHI(I,J,K,1:5))

               A = CMPLX(0.0D0,0.0D0,KIND=DBL)
               W = 0.0D0
            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO
  END SELECT
END SUBROUTINE SE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                     Calculates the discrete forward Fourier transform                            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE FFT()
  USE FFTW3, ONLY : FFTW_EXECUTE_DFT
  USE CGPE_DATA, ONLY : PHIF, PHI, SPIN
  USE FFTW_DATA
 
  IMPLICIT NONE
 
  INTEGER :: L, MS

  MS = 2*SPIN+1
  DO L = 1, MS
     FFTFXYZ(:,:,:) = PHI(:,:,:,L)
     CALL FFTW_EXECUTE_DFT(PLANFXYZ,FFTFXYZ,FFTBXYZ)
     PHIF(:,:,:,L) = FFTBXYZ(:,:,:)
  END DO
END SUBROUTINE FFT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         Calculates the discrete backward Fourier transform of component wavefunctions            !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE BFT()
  USE BASIC_DATA, ONLY : NX, NY, NZ
  USE CGPE_DATA, ONLY : PHIF, PHI, SPIN
  USE FFTW3, ONLY : FFTW_EXECUTE_DFT
  USE FFTW_DATA

  IMPLICIT NONE

  INTEGER :: L, MS

  MS = 2*SPIN+1
  DO L = 1, MS
     FFTBXYZ(:,:,:) = PHIF(:,:,:,L)
     CALL FFTW_EXECUTE_DFT(PLANBXYZ,FFTBXYZ,FFTFXYZ)
     PHI(:,:,:,L) = FFTFXYZ(:,:,:)/DBLE(NX*NY*NZ)
  END DO
END SUBROUTINE BFT

