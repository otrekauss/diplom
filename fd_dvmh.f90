!===============================================================================
! Модуль конечно-разностных схем для DVMH
! Заменяет ручное управление CUDA ядрами на параллельные циклы DVM
!===============================================================================
MODULE FD_DVMH
  USE CONSTANTS
  USE VAR_DVMH
  IMPLICIT NONE
  
CONTAINS

  !------------------------------------------------------
  SUBROUTINE SIMPLER
    USE MD_DVMH, ONLY: MHD, ComputeRotB
    USE SOLVER_DVMH, ONLY: AdiSolver
    USE USER_DVMH, ONLY: ComputeSourceU1, ComputeSourceU2, ComputeSourceU3, &
                         ComputeSourceP, ComputeSourcePC, ComputeSourceT 
    IMPLICIT NONE
    INTEGER :: NUM_ITERS, OUT_ITER
    
    SMAX = 1.0_8
    CALL ComputeRotB
    
    !$DVM TASK SHARED
    CALL ComputeFOnBot(U1)
    CALL ComputeFOnBot(U3)
    CALL ComputeFOnBot(T)
    !$DVM END TASK
    
    DO WHILE (SMAX > TOL_SMAX)
      
      !-----------------Step 1------------------
      CALL ComputeSourceU1
      CALL CofU1
      
      CALL ComputeSourceU2
      CALL CofU2
      
      CALL ComputeSourceU3
      CALL CofU3
      
      CALL ComputeSourceP
      CALL CofP
      NUM_ITERS = ADISolver(P, N1_core, N1-1, 2, 2, NTIMES_P)
      
      !-----------------Step 2---------------
      U1tmp = U1
      U2tmp = U2
      U3tmp = U3
      
      CALL ComputeSourceU1
      CALL CofU1
      NUM_ITERS = AdiSolver(U1tmp, N1_core + 1, N1-1, 2, 2, NTIMES_UVW)
      CALL ComputeFOnBot(U1tmp)
      
      CALL ComputeSourceU2
      CALL CofU2
      NUM_ITERS = ADISolver(U2tmp, N1_core, N1-1, 3, 2, NTIMES_UVW)
      
      CALL ComputeSourceU3
      CALL CofU3
      NUM_ITERS = ADISolver(U3tmp, N1_core, N1-1, 2, 2, NTIMES_UVW)
      CALL ComputeFOnBot(U3tmp)
      
      !-----------------Step 3---------------
      CALL ComputeSourcePC
      CALL CofPC
      NUM_ITERS = ADISolver(PC, N1_core, N1-1, 2, 2, NTIMES_PC)
      
      !-----------------Step 4---------------
      U1 = U1tmp
      U2 = U2tmp
      U3 = U3tmp
      CALL CorrectU
      
      CALL ComputeSmax
      
    END DO 
    
    !-----------------Step 5---------------
    CALL MHD
    CALL ComputeRotB
    
    CALL ComputeSourceT
    CALL CofT
    NUM_ITERS = ADISolver(T, N1_core, N1-1, 2, 2, NTIMES_T)
    CALL ComputeFOnBot(T)
    
  END SUBROUTINE SIMPLER
  
  !------------------------------------------------------
  SUBROUTINE ComputeFOnBot(F)
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:,:), INTENT(INOUT) :: F
    INTEGER :: I, J
    
    !$DVM TASK SHARED PRIVATE(I,J)
    DO J = 1, N2
      DO I = 1, N1_core
        F(I,J,1) = 0.0_8
      ENDDO
    ENDDO
    !$DVM END TASK
    
  END SUBROUTINE ComputeFOnBot
  
  !------------------------------------------------------
  SUBROUTINE CofU1
    USE MD_DVMH, ONLY: CofAimAipU1, CofAjmAjpU1, CofAkmAkpU1
    IMPLICIT NONE
    
    !$DVM TASK SHARED
    CALL CofAimAipU1
    CALL CofAjmAjpU1
    CALL CofAkmAkpU1
    CALL CofConApU1hatDU1
    CALL CofConU1
    !$DVM END TASK
    
  END SUBROUTINE CofU1
  
  !------------------------------------------------------
  SUBROUTINE CofU2
    USE MD_DVMH, ONLY: CofAimAipU2, CofAjmAjpU2, CofAkmAkpU2
    IMPLICIT NONE
    
    !$DVM TASK SHARED
    CALL CofAimAipU2
    CALL CofAjmAjpU2
    CALL CofAkmAkpU2
    CALL CofConApU2hatDU2
    CALL CofConU2
    !$DVM END TASK
    
  END SUBROUTINE CofU2
  
  !------------------------------------------------------
  SUBROUTINE CofU3
    IMPLICIT NONE
    
    !$DVM TASK SHARED
    CALL CofU3_core
    CALL FF1N1(DU3, N1)
    CALL FF1N1(U3hat, N1)
    !$DVM END TASK
    
  END SUBROUTINE CofU3
  
  !------------------------------------------------------
  SUBROUTINE CofU3_core
    IMPLICIT NONE
    INTEGER :: I, J, K
    
    !$DVM TASK SHARED PRIVATE(I,J,K)
    DO K = 1, N3
      DO J = 1, N2
        DO I = N1_core, N1-1
          ! Вычисления для U3
          DU3(I,J,K) = 0.0_8  ! Заглушка - заменить на реальные вычисления
        ENDDO
      ENDDO
    ENDDO
    !$DVM END TASK
    
  END SUBROUTINE CofU3_core
  
  !------------------------------------------------------
  SUBROUTINE CofP
    IMPLICIT NONE
    ! Коэффициенты для уравнения давления
    !$DVM TASK SHARED
    ! Вычисления...
    !$DVM END TASK
  END SUBROUTINE CofP
  
  !------------------------------------------------------
  SUBROUTINE CofPC
    IMPLICIT NONE
    ! Коэффициенты для корректировки давления
    !$DVM TASK SHARED
    ! Вычисления...
    !$DVM END TASK
  END SUBROUTINE CofPC
  
  !------------------------------------------------------
  SUBROUTINE CofT
    IMPLICIT NONE
    ! Коэффициенты для уравнения температуры
    !$DVM TASK SHARED
    ! Вычисления...
    !$DVM END TASK
  END SUBROUTINE CofT
  
  !------------------------------------------------------
  SUBROUTINE CorrectU
    IMPLICIT NONE
    ! Коррекция скоростей
    !$DVM TASK SHARED
    ! Вычисления...
    !$DVM END TASK
  END SUBROUTINE CorrectU
  
  !------------------------------------------------------
  SUBROUTINE ComputeSmax
    IMPLICIT NONE
    REAL(8) :: SMAX_LOCAL
    INTEGER :: I, J, K
    
    SMAX_LOCAL = 0.0_8
    
    !$DVM TASK SHARED PRIVATE(I,J,K) REDUCTION(MAX:SMAX_LOCAL)
    DO K = 1, N3
      DO J = 1, N2
        DO I = 1, N1-1
          SMAX_LOCAL = MAX(SMAX_LOCAL, ABS(DU1(I,J,K)), ABS(DU2(I,J,K)), ABS(DU3(I,J,K)))
        ENDDO
      ENDDO
    ENDDO
    !$DVM END TASK
    
    SMAX = SMAX_LOCAL
    
  END SUBROUTINE ComputeSmax
  
END MODULE FD_DVMH
